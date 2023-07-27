#!/usr/bin/env python

import sys
from igraph import *

def read_graph_gml(gml_fname):
    G = read(gml_fname)
    G.vs["id"] = [int(v["id"]) for v in G.vs]
    G.vs["readlen"] = [int(v["readlen"]) for v in G.vs]
    G.es["direction"] = [int(e["direction"]) for e in G.es]
    G.es["suffix"] = [int(e["suffix"]) for e in G.es]
    G.es["prefix"] = [int(e["prefix"]) for e in G.es]
    return G

def identify_tips(G):
    roots = G.vs.select(_indegree_eq=1)
    branches = G.vs.select(_indegree_ge=3)
    tips = G.es.select(_between=(roots, branches))
    return tips

def main(input_gml, output_gml):
    G = read_graph_gml(input_gml)
    tips = identify_tips(G)
    tips.delete()
    G.write_gml(output_gml)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("usage: {} <graph.gml> <graph.notips.gml>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    input_gml = sys.argv[1]
    output_gml = sys.argv[2]

    main(input_gml, output_gml)

    sys.exit(0)
