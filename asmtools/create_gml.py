#!/usr/bin/env python

import sys
from igraph import *

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

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

def main(elba_paf_fname, reads_fname, elba_gml_fname):

    readseqs, readnames, namemap = read_fasta(reads_fname)
    readlens = [len(s) for s in readseqs]

    overlaps, attrs = read_paf(elba_paf_fname, readlens, namemap)
    S = create_overlap_graph(readlens, readnames, overlaps, attrs)

    S.write_gml(elba_gml_fname)

if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <elba.paf> <reads.fa> <elba.gml>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)
    elba_paf_fname = sys.argv[1]
    reads_fname = sys.argv[2]
    elba_gml_fname = sys.argv[3]

    main(elba_paf_fname, reads_fname, elba_gml_fname)

    sys.exit(0)
