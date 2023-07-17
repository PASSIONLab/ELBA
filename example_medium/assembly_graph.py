#!/usr/bin/env python

import sys
import getopt
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

neighborhood = 3
output_prefix = "asm"

def usage():
    global neighborhood, output_prefix
    sys.stderr.write("Usage: {} [options] <elba.string.paf> <elba.overlap.paf> <reads.fa>\n".format(sys.argv[0]))
    sys.stderr.write("Options: -l INT   branch/bridge neighborhood value [{}]\n".format(neighborhood))
    sys.stderr.write("         -o STR   output file prefix [{}]\n".format(output_prefix))
    sys.stderr.write("         -h       help message\n")
    return -1

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

def read_paf(paf_fname, readlens, namemap):
    overlaps = []
    attrs = []
    types = (str, int, int, int, lambda x: int(x[0] == '-'), str, int, int, int)
    for line in open(paf_fname, "r"):
        tokens = line.rstrip().split()[:9]
        nameQ, lenQ, begQ, endQ, strand, nameT, lenT, begT, endT = [t(tok) for t, tok in zip(types, tokens)]
        assert nameQ in namemap and nameT in namemap
        idQ, idT = namemap[nameQ], namemap[nameT]
        assert lenQ == readlens[idQ] and lenT == readlens[idT]
        overlaps.append((idQ, idT))
        attrs.append((begQ, endQ, strand, begT, endT))
    return overlaps, attrs

def create_overlap_graph(readlens, readnames, overlaps, attrs):
    n = len(readlens)
    G = Graph(n, directed=True)
    G.vs["readlen"] = readlens
    G.vs["readname"] = readnames
    edges, directions, suffixes, prefixes = [], [], [], []
    for overlap, attr in zip(overlaps, attrs):
        idQ, idT = overlap
        if idQ >= idT: continue
        begQ, endQ, strand, begT, endT = attr
        lenQ, lenT = readlens[idQ], readlens[idT]
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

def assemble_chain(chain, readseqs):
    parts = []
    for readid, pre, strand in chain:
        s = readseqs[readid]
        if strand == 1: s = revcomp(s)
        parts.append(s[:pre])
    return "".join(parts)

def identify_tips(G):
    roots = G.vs.select(_indegree_eq=1)
    branches = G.vs.select(_indegree_ge=3)
    tips = G.es.select(_between=(roots, branches))
    return tips

def identify_bridges(G):
    branches = G.vs.select(_indegree_ge=3)
    bridges = []
    targets = set()
    for branch in branches:
        for t in branch.successors():
            if t.index in targets:
                bridges.append(t.index)
            else:
                targets.add(t.index)
    return bridges

def extend_chain(G, root, visited, l):
    chain = [root]
    cur = root
    for i in range(l-1):
        front = G.vs[cur].successors()
        neighs = [x.index for x in front if x.index not in visited]
        if len(neighs) != 1: break
        visited.add(cur)
        cur = neighs.pop()
        chain.append(cur)
    return chain


def extend_chain(G, root, visited, l):
    chain = [root]
    cur = root
    for i in range(l-1):
        front = G.vs[cur].successors()
        neighs = [x.index for x in front if x.index not in visited]
        if len(neighs) != 1: break
        visited.add(cur)
        cur = neighs.pop()
        chain.append(cur)
    return chain

def bridge_neighborhood_chains(G, bridge, l):

    b = G.vs[bridge]

    if b.indegree() != 2: return None

    westend, eastend = [x.index for x in b.successors()]

    if G.vs[westend].indegree() != 3 or G.vs[eastend].indegree() != 3: return None

    visited = {bridge, westend, eastend}

    west_chain_roots = [x.index for x in G.vs[westend].successors() if x.index not in visited]
    east_chain_roots = [x.index for x in G.vs[eastend].successors() if x.index not in visited]

    if len(west_chain_roots) != 2 or len(east_chain_roots) != 2: return None

    northwest_root, southwest_root = west_chain_roots
    northeast_root, southeast_root = east_chain_roots

    northwest_chain = extend_chain(G, northwest_root, visited, l)
    southwest_chain = extend_chain(G, southwest_root, visited, l)
    northeast_chain = extend_chain(G, northeast_root, visited, l)
    southeast_chain = extend_chain(G, southeast_root, visited, l)

    chains = [northwest_chain, southwest_chain, northeast_chain, southeast_chain]

    for chain in chains:
        if len(chain) != l: return None

    return chains

def write_contigs(S, outprefix, readseqs):
    chains = generate_contig_chains(S.copy())
    contigs = [assemble_chain(chain, readseqs) for chain in chains]
    with open("{}.contigs.fa".format(outprefix), "w") as f:
        for i, contig in enumerate(contigs):
            f.write(">contig{}\n{}\n".format(i+1, contig))

def main(argc, argv):
    global neighborhood, output_prefix

    try: opts, args = getopt.gnu_getopt(argv[1:], "l:o:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-l": neighborhood = int(a)
        elif o == "-o": output_prefix = a
        elif o == "-h": return usage()


    if len(args) != 3:
        return usage()

    string_paf_fname = args[0]
    overlap_paf_fname = args[1]
    fasta_fname = args[2]

    readseqs, readnames, namemap = read_fasta(fasta_fname)
    readlens = [len(s) for s in readseqs]

    overlaps, attrs = read_paf(string_paf_fname, readlens, namemap)
    S = create_overlap_graph(readlens, readnames, overlaps, attrs)
    S.write_gml("{}.S.gml".format(output_prefix))
    write_contigs(S, output_prefix, readseqs)

    overlaps, attrs = read_paf(overlap_paf_fname, readlens, namemap)
    R = create_overlap_graph(readlens, readnames, overlaps, attrs)
    R.write_gml("{}.R.gml".format(output_prefix))

    itr = 1

    while True:

        # step one: remove tips
        tips = identify_tips(S)
        tips_found = len(tips)
        tips.delete()
        sys.stderr.write("Iteration {}: Found and removed {} tips\n".format(itr, tips_found//2))

        # step two: identify bridges
        bridges = identify_bridges(S)
        sys.stderr.write("Iteration {}: Found {} bridges\n".format(itr, len(bridges)))

        # step 3: identify edges between bridges and delete them
        double_bridges = S.es.select(_within=bridges)

        if len(double_bridges) == 0: break
        else:
            sys.stderr.write("Iteration {}: Found and removed {} double bridges\n".format(itr, len(double_bridges)))
            double_bridges.delete()

        itr += 1

    while True:

        S.vs["bridge"] = 0
        bridges = identify_bridges(S)
        S.vs[bridges]["bridge"] = 1
        S.es["added"] = 0
        S.es["remove"] = 0
        for bridge in bridges:
            bridge_chains = bridge_neighborhood_chains(S, bridge, neighborhood)
            if bridge_chains is None: continue
            counts = []
            for i in range(3):
                for j in range(i+1, 4):
                    num = len(R.es.select(_between=(bridge_chains[i], bridge_chains[j])))//2
                    counts.append(num)
            if counts[0] > 0 and counts[5] > 0 and sum(counts[1:5]) == 0:
                toremove = S.es.select(_from=bridge)
                for e in toremove: e["remove"] = 1
                toremove = S.es.select(_to=bridge)
                for e in toremove: e["remove"] = 1

        es = S.es.select(remove_eq=1)
        if len(es) == 0: break
        es.delete()

        tips = identify_tips(S)
        tips.delete()

    S.write_gml("{}.S.cleaned.gml".format(output_prefix))
    write_contigs(S, "{}.cleaned".format(output_prefix), readseqs)

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
