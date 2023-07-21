#!/usr/bin/env python

import sys
import getopt
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

output_prefix = "asm"

def usage():
    global output_prefix
    sys.stderr.write("Usage: {} [options] <elba.string.gml> <reads.fa>\n".format(sys.argv[0]))
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

def write_contigs(S, outprefix, readseqs):
    chains = generate_contig_chains(S.copy())
    contigs = [assemble_chain(chain, readseqs) for chain in chains]
    with open("{}.contigs.fa".format(outprefix), "w") as f:
        for i, contig in enumerate(contigs):
            f.write(">contig{}\n{}\n".format(i+1, contig))

def read_graph_gml(gml_fname):
    G = read(gml_fname)
    G.vs["id"] = [int(v["id"]) for v in G.vs]
    G.vs["readlen"] = [int(v["readlen"]) for v in G.vs]
    G.es["direction"] = [int(e["direction"]) for e in G.es]
    G.es["suffix"] = [int(e["suffix"]) for e in G.es]
    G.es["prefix"] = [int(e["prefix"]) for e in G.es]
    return G

def main(argc, argv):
    global output_prefix

    try: opts, args = getopt.gnu_getopt(argv[1:], "o:h")
    except getopt.GetoptError as err:
        sys.stderr.write("error: {}\n".format(err))
        return usage()

    for o, a in opts:
        if o == "-o": output_preix = a
        elif o == "-h": return usage()

    if len(args) != 2:
        return usage()

    string_gml_fname = args[0]
    fasta_fname = args[1]

    readseqs, readnames, namemap = read_fasta(fasta_fname)
    readlens = [len(s) for s in readseqs]

    S = read_graph_gml(string_gml_fname)

    write_contigs(S, output_prefix, readseqs)

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
