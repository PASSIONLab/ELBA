#!/usr/bin/env python

import sys
import subprocess as sp
import shutil
from pathlib import Path
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def revcomp(s): return s.translate(comp_tab)[::-1]

def read_fasta(fasta_fname):
    seqs, names, seq, name = [], [], [], ""
    for line in open(fasta_fname):
        if line[0] == ">":
            if len(seq) > 0:
                seqs.append("".join(seq))
                names.append(name)
            name = line.lstrip(">").split()[0].rstrip()
            seq = []
        else: seq.append(line.rstrip())
    if len(seq) > 0:
        seqs.append("".join(seq))
        names.append(name)
    return seqs, names

def create_name_map(names):
    return {names[i] : i for i in range(len(names))}

def create_elba_graph(paf_fname, seqs, names):
    tuples = []
    n = len(seqs)
    assert n == len(names)
    namemap = create_name_map(names)
    types = (str, int, int, int, lambda x: int(x[0] == '-'), str, int, int, int)
    for line in open(paf_fname, "r"):
        tokens = line.rstrip().split()[:9]
        nameQ, lenQ, begQ, endQ, strand, nameT, lenT, begT, endT = [t(tok) for t, tok in zip(types, tokens)]
        assert nameQ in namemap and nameT in namemap
        idQ, idT = namemap[nameQ], namemap[nameT]
        if idQ >= idT: continue
        assert lenQ == len(seqs[idQ]) and lenT == len(seqs[idT])
        tuples.append((idQ, lenQ, begQ, endQ, strand, idT, lenT, begT, endT))
    G = Graph(n, directed=True)
    G.vs["readlen"] = [len(s) for s in seqs]
    G.vs["readname"] = names
    edges, directions, suffixes, prefixes = [], [], [], []
    for idQ, lenQ, begQ, endQ, strand, idT, lenT, begT, endT in tuples:
        begTr = begT if strand == 0 else lenT - endT
        endTr = endT if strand == 0 else lenT - begT
        if begQ <= begTr and lenQ - endQ <= lenT - endTr:
            pass
        elif begQ >= begTr and lenQ - endQ >= lenT - endTr:
            pass
        else:
            edges += [(idQ, idT), (idT, idQ)]
            if begQ > begTr:
                lengths = [(lenT - endTr) - (lenQ - endQ), begQ - begTr]
                suffixes += lengths
                prefixes += lengths[::-1]
                if strand == 0:
                    directions += [1, 2]
                else:
                    directions += [0, 0]
            else:
                lengths = [begTr - begQ, (lenQ - endQ) - (lenT - endTr)]
                suffixes += lengths
                prefixes += lengths[::-1]
                if strand == 0:
                    directions += [2, 1]
                else:
                    directions += [3, 3]
    G.add_edges(edges)
    G.es["direction"] = directions
    G.es["suffix"] = suffixes
    G.es["prefix"] = prefixes
    return G

def assemble_chain(chain, readseqs):
    parts = []
    info = []
    for readid, pre, strand in chain:
        s = readseqs[readid]
        if strand == 1: s = revcomp(s)
        parts.append(s[:pre])
        a = 0 if strand == 0 else len(s)-1
        b = pre if strand == 0 else len(s)-pre
        info.append((readid, a, b))
    return "".join(parts), info

def generate_contig_chains(S):
    G = S.copy()
    branches = G.vs.select(_indegree_ge=3)
    branch_edges = G.es.select(_between=(branches, G.vs))
    G.delete_edges(branch_edges)
    contig_components = G.components()
    contig_chains = []
    for contig_component in contig_components:
        q = len(contig_component)
        if q <= 1: continue
        vertices = G.vs.select(contig_component)
        roots = vertices.select(_indegree_eq=1)
        if len(roots) != 2: continue
        assert len(roots) == 2
        s, t = roots
        topsort = G.bfs(s)[0]
        if topsort[0] != s.index or topsort[-1] != t.index: continue
        chain = []
        last_dir = None
        for i in range(q-1):
            e = G.es[G.get_eid(topsort[i], topsort[i+1])]
            strand = int((e["direction"]>>1)&1)
            chain.append((e.source, e["prefix"], strand))
            last_dir = e["direction"]
        chain.append((t.index, t["readlen"], 1 - int(last_dir&1)))
        contig_chains.append(chain)
    return contig_chains

def read_vector(fname):
    vec = []
    for line in open(fname, "r"):
        vec.append(int(line.rstrip().split()[1]))
    return vec

def write_neighbors(fname, vertices, G):
    with open(fname, "w") as f:
        for i in sorted(vertices):
            u = G.vs[i]
            ns = sorted([v.index for v in u.successors()])
            for v in ns: f.write("{}\t{}\n".format(u.index, v))

def get_removed_edges(before, after):
    after_edges = {}
    removed_edges = {}
    for e in after.es:
        u, v = e.tuple
        if not u in after_edges: after_edges[u] = set()
        if not v in after_edges[u]: after_edges[u].add(v)
    for e in before.es:
        u, v = e.tuple
        if not u in after_edges or not v in after_edges[u]:
            if not u in removed_edges: removed_edges[u] = set()
            if not v in removed_edges[u]: removed_edges[u].add(v)
    return removed_edges

def expand_edges(edges):
    expanded = []
    for u in sorted(edges.keys()):
        vs = sorted(list(edges[u]))
        for v in vs:
            expanded.append((u, v))
    return expanded

def get_bridges(G):
    visited = set()
    bridges = set()
    for triple in G.vs.select(_indegree_ge=3):
        for u in triple.successors():
            if not u.index in visited:
                visited.add(u.index)
            else:
                bridges.add(u.index)
    return bridges

def get_experiments(params_fname):
    experiments = {}
    with open(params_fname, "r") as f:
        header = next(f).rstrip()
        assert header == "L,U,k,A,B,f,s,x,c,reads_fname,ref_fname,outdir"
        for line in f.readlines():
            L,U,k,A,B,fr,s,x,c,fasta_fname,ref_fname,outdir = (t(tok) for t,tok in zip((int,int,int,int,int,float,int,int,float,str,str,str), line.rstrip().split(",")))
            fasta_path = Path(fasta_fname).resolve()
            ref_path = Path(ref_fname).resolve()
            outdir_path = Path(outdir).resolve()
            assert fasta_path.is_file()
            assert ref_path.is_file()
            assert outdir_path.is_dir()
            if not L in experiments: experiments[L] = {}
            if not U in experiments[L]: experiments[L][U] = {}
            if not k in experiments[L][U]: experiments[L][U][k] = []
            experiments[L][U][k].append((A, B, fr, s, x, c, fasta_path, ref_path, outdir_path))
    return experiments

def process_exp(exp, seqs, names):
    fasta_path, ref_path, outdir_path = exp[-3:]
    pathdict = {}
    for path in outdir_path.iterdir():
        pathdict[path.name] = str(path.resolve())
    for name in ["A.paf", "B.paf", "C.paf", "elba.overlap.paf", "elba.string.paf", "contained_reads.txt", "bad_reads.txt"]:
        assert name in pathdict
    A = create_elba_graph(pathdict["A.paf"], seqs, names)
    B = create_elba_graph(pathdict["B.paf"], seqs, names)
    R = create_elba_graph(pathdict["elba.overlap.paf"], seqs, names)
    C = create_elba_graph(pathdict["C.paf"], seqs, names)
    S = create_elba_graph(pathdict["elba.string.paf"], seqs, names)
    bad_reads = sorted(read_vector(pathdict["bad_reads.txt"]))
    contained_reads = sorted(read_vector(pathdict["contained_reads.txt"]))
    bad_edges = expand_edges(get_removed_edges(A, B))
    contained_edges = expand_edges(get_removed_edges(R, C))
    transitive_edges = expand_edges(get_removed_edges(C, S))
    branching_reads = sorted(S.vs.select(_indegree_ge=3).indices)
    bridge_reads = sorted(list(get_bridges(S)))

    chains = generate_contig_chains(S)
    contigs = []

    for chain in chains:
        contig, info = assemble_chain(chain, seqs)
        contigs.append((contig, info))

    contigpath = outdir_path.joinpath("exp.contigs.fa").resolve()
    assert not contigpath.exists()

    with open(str(contigpath), "w") as f:
        for i, (contig, info) in enumerate(contigs):
            f.write(">contig{}\n{}\n".format(i+1, contig))

    quastpath = outdir_path.joinpath("quastdir").resolve()
    assert not quastpath.exists()
    quast_cmd = ["quast.py", "-o", str(quastpath), "-r", str(ref_path), str(contigpath)]
    proc = sp.Popen(quast_cmd, stdout=sp.PIPE)
    proc.wait()

    assert quastpath.is_dir()
    reportpath = quastpath.joinpath("transposed_report.tsv").resolve()
    assert reportpath.is_file()

    header = ["L", "U", "k", "A", "B", "f", "s", "x", "c"]
    values = [str(item) for item in exp[:9]]

    with open(str(reportpath), "r") as f:
        header += next(f).rstrip().split("\t")
        values += next(f).rstrip().split("\t")

    statspath = outdir_path.joinpath("elba.asm.stats.tsv").resolve()
    assert not statspath.exists()

    with open(str(statspath), "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(values) + "\n")

    procpath = outdir_path.joinpath("elba.graph.stats.txt").resolve()
    assert not procpath.exists()

    with open(str(procpath), "w") as f:
        f.write("# ... - {} reads total\n".format(len(S.vs)))
        f.write("# ^AA - {} bad reads\n".format(len(bad_reads)))
        f.write("# ^BB - {} contained reads\n".format(len(contained_reads)))
        f.write("# ^CC - {} branching reads\n".format(len(branching_reads)))
        f.write("# ^DD - {} bridge reads\n".format(len(bridge_reads)))
        f.write("# ^EE - {} bad overlaps\n".format(len(bad_edges)//2))
        f.write("# ^FF - {} contained overlaps\n".format(len(contained_edges)//2))
        f.write("# ^GG - {} transitive overlaps\n".format(len(transitive_edges)//2))
        f.write("# ^HH - contig info\n")
        f.write("# ^II - read assignments\n")

        f.write("AA\t{} bad reads\n".format(len(bad_reads)))
        for i in bad_reads: f.write("AA\t{}\n".format(i))

        f.write("BB\t{} contained reads\n".format(len(contained_reads)))
        for i in contained_reads: f.write("BB\t{}\n".format(i))

        f.write("CC\t{} branching reads\n".format(len(branching_reads)))
        for i in branching_reads: f.write("CC\t{}\n".format(i))

        f.write("DD\t{} bridge reads\n".format(len(bridge_reads)))
        for i in bridge_reads: f.write("DD\t{}\n".format(i))

        f.write("EE\t{} bad overlaps\n".format(len(bad_edges)//2))
        for u,v in bad_edges: f.write("EE\t{}\t{}\n".format(u,v))

        f.write("FF\t{} contained overlaps\n".format(len(contained_edges)//2))
        for u,v in contained_edges: f.write("FF\t{}\t{}\n".format(u,v))

        f.write("GG\t{} transitive overlaps\n".format(len(transitive_edges)//2))
        for u,v in transitive_edges: f.write("GG\t{}\t{}\n".format(u,v))

        f.write("HH\t{contig-id}\t{read-id}\t{read-start}\t{read-end}\t{contig-start}\t{contig-end}\n")

        read_assignments = {}

        for i, (contig, info) in enumerate(contigs):
            contigpos = 0
            for readid, a, b in info:
                read_assignments[readid] = i
                f.write("HH\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i+1, readid, a, b, contigpos, contigpos + abs(a-b)))
                contigpos += abs(a-b)

        f.write("II\t{read-id}\t{contig-id}\n")
        for i in range(len(S.vs)):
            if not i in read_assignments:
                f.write("II\t{}\t(null)\n".format(i))
            else:
                f.write("II\t{}\t{}\n".format(i, read_assignments[i]+1))

    spath = outdir_path.joinpath("string_graph.gml").resolve()
    assert not spath.exists()

    S.vs["contigid"] = 0
    for readid, contigid in read_assignments.items():
        S.vs[readid]["contigid"] = contigid+1

    S.write_gml(str(spath))
    return header, values

def main(argc, argv):
    experiments = get_experiments(argv[1])
    vseqs, vnames = {}, {}
    valuelist = []
    for L in experiments:
        for U in experiments[L]:
            for k in experiments[L][U]:
                for exp in experiments[L][U][k]:
                    fasta_path, ref_path, outdir_path = exp[-3:]
                    if not fasta_path.name in vseqs:
                        seqs, names = read_fasta(str(fasta_path))
                        vseqs[fasta_path.name] = seqs
                        vnames[fasta_path.name] = names
                    explist = [L, U, k] + list(exp)
                    header, values = process_exp(explist, vseqs[fasta_path.name], vnames[fasta_path.name])
                    valuelist.append(values)
    with open(Path.cwd().joinpath("elba.asm.stats.csv"), "w") as f:
        f.write(",".join(header) + "\n")
        for values in valuelist:
            f.write(",".join(values) + "\n")

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
