#!/usr/bin/env python

import sys
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def revcomp(s): return s.translate(comp_tab)[::-1]

def read_fasta(fasta_fname):
    readseqs = list()
    readnames = list()
    seq = []
    readname = ""
    for line in open(fasta_fname, "r"):
        if line.startswith(">"):
            if len(seq) > 0:
                readseqs.append("".join(seq))
                readnames.append(readname)
            readname = line.lstrip(">").split()[0].rstrip()
            seq = []
        else:
            seq.append(line.rstrip())
    if len(seq) > 0:
        readseqs.append("".join(seq))
        readnames.append(readname)
    namemap = dict()
    for i, readname in enumerate(readnames):
        namemap[readname] = i
    return readseqs, readnames, namemap

def read_graph_gml(gml_fname):
    G = read(gml_fname)
    G.vs["id"] = [int(v["id"]) for v in G.vs]
    G.vs["readlen"] = [int(v["readlen"]) for v in G.vs]
    G.es["direction"] = [int(e["direction"]) for e in G.es]
    G.es["suffix"] = [int(e["suffix"]) for e in G.es]
    G.es["prefix"] = [int(e["prefix"]) for e in G.es]
    return G

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
        assert len(roots) == 2
        s, t = roots
        topsort = G.bfs(s)[0]
        assert topsort[0] == s.index and topsort[-1] == t.index
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

def get_mapped_regions(e, readseqs):
    qid, tid = e.tuple
    qs, ts = readseqs[qid], readseqs[tid]
    qstrand = int((e["direction"]>>1)&1)
    tstrand = 1 - int(e["direction"]&1)
    if qstrand == 1: qs = revcomp(qs)
    if tstrand == 1: ts = revcomp(ts)
    qregion = qs[e["prefix"]:]
    tregion = ts[:len(ts)-e["suffix"]]
    return qregion, tregion

def similarity(s1, s2, k):
    s1kmers = set(hash(s1[i:i+k]) for i in range(len(s1)-k+1))
    s2kmers = set(hash(s2[i:i+k]) for i in range(len(s2)-k+1))
    shared = len(s1kmers.intersection(s2kmers))
    total = len(s1kmers.union(s2kmers))
    return shared / total

def assemble_chain(chain, readseqs):
    parts = []
    for readid, pre, strand in chain:
        s = readseqs[readid]
        if strand == 1: s = revcomp(s)
        parts.append(s[:pre])
    return "".join(parts)

def write_contigs(S, fname, readseqs):
    chains = generate_contig_chains(S.copy())
    H = S.copy()
    H.vs["contigid"] = 0
    with open(fname, "w") as f:
        for i, chain in enumerate(chains):
            contig = assemble_chain(chain, readseqs)
            vertices = [item[0] for item in chain]
            H.vs[vertices]["contigid"] = i+1
            f.write(">contig{}\t{}-{}\n{}\n".format(i+1, vertices[0], vertices[-1], contig))
    H.write_gml("H.gml")

def main(elba_gml_fname, reads_fname, elba_contigs_fname):

    readseqs, readnames, namemap = read_fasta(reads_fname)
    readlens = [len(s) for s in readseqs]

    S = read_graph_gml(elba_gml_fname)
    write_contigs(S, elba_contigs_fname, readseqs)


if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <elba.gml> <reads.fa> <elba.contigs.fa>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    elba_gml_fname = sys.argv[1]
    reads_fname = sys.argv[2]
    elba_contigs_fname = sys.argv[3]

    main(elba_gml_fname, reads_fname, elba_contigs_fname)

    sys.exit(0)
