#!/usr/bin/env python

import sys
from igraph import Graph

def read_fasta(fasta):
    read_lengths = []
    read_names = []
    with open(fasta, "r") as f:
        seqlen = 0
        name = ""
        for line in f.readlines():
            if line.startswith(">"):
                if seqlen > 0:
                    read_lengths.append(seqlen)
                    read_names.append(name)
                    seqlen = 0
                name = line.lstrip(">").rstrip()
            else:
                seqlen += len(line.rstrip())
        if seqlen > 0:
            read_lengths.append(seqlen)
            read_names.append(name)
    assert len(read_lengths) == len(read_names)
    return read_names, read_lengths

def get_vertices(fasta_fname):
    vertices = []
    read_names, read_lengths = read_fasta(fasta_fname)
    read_name_map = {}
    for i in range(len(read_names)):
        read_name = read_names[i]
        read_length = read_lengths[i]
        names = read_name.split(None, 1)
        if len(names) == 2:
            name, comment = names
        else:
            name = names[0]
            comment = ""
        read_name_map[name] = i
        vertices.append({"name" : name, "comment" : comment, "length" : read_length})
    return vertices, read_name_map, read_lengths

def get_edges(paf_fname, vertices, read_name_map, read_lengths):
    edges = []
    for line in open(paf_fname, "r"):
        items = line.rstrip().split()
        source_name, target_name = items[0], items[5]
        source_length, target_length = int(items[1]), int(items[6])
        assert source_length == read_lengths[read_name_map[source_name]]
        assert target_length == read_lengths[read_name_map[target_name]]
        source_beg, target_beg = int(items[2]), int(items[7])
        source_end, target_end = int(items[3]), int(items[8])
        strand = items[4]
        edges.append({"source" : source_name, "target" : target_name, "source_beg" : source_beg, "source_end" : source_end, "target_beg" : target_beg, "target_end" : target_end, "strand" : strand})
    return edges

def main(fasta_fname, paf_fname, gml_fname):

    vertices, read_name_map, read_lengths = get_vertices(fasta_fname)

    edges = get_edges(paf_fname, vertices, read_name_map, read_lengths)

    graph = Graph.DictList(vertices, edges, directed=True, vertex_name_attr="name", edge_foreign_keys=('source', 'target'), iterative=False)

    del graph.es["source"]
    del graph.es["target"]

    f = gml_fname if gml_fname != "stdout" else sys.stdout

    graph.write_gml(f, creator="Gabriel Raulet; paf2gml.py")

    return 0

if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <reads.fa> <ava.paf> <ava.gml | stdout>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    fasta_fname = sys.argv[1]
    paf_fname = sys.argv[2]
    gml_fname = sys.argv[3]

    retval = main(fasta_fname, paf_fname, gml_fname)

    sys.exit(retval)
