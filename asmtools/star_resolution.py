#!/usr/bin/env python

import sys
from igraph import *
from collections import Counter

def read_graph_gml(gml_fname):
    G = read(gml_fname)
    G.vs["id"] = [int(v["id"]) for v in G.vs]
    G.vs["readlen"] = [int(v["readlen"]) for v in G.vs]
    G.es["direction"] = [int(e["direction"]) for e in G.es]
    G.es["suffix"] = [int(e["suffix"]) for e in G.es]
    G.es["prefix"] = [int(e["prefix"]) for e in G.es]
    return G

def star_resolution(S, R):
    stars = []
    for u in S.vs.select(_indegree_eq=3):
        isstar = True
        for v in u.successors():
            if v.indegree() != 2:
                isstar = False
        if isstar:
            stars.append(u.index)

    S.vs["star"] = 0
    S.vs[stars]["star"] = 1

    star_edges = []
    star_verts = []

    for u in S.vs.select(star=1):
        neighs = set([a.index for a in u.successors()])
        star_arcs = R.es.select(_within=neighs)
        if len(star_arcs) == 2:
            for e in R.es.select(_within=neighs):
                star_edges.append(e.tuple)
            arcverts = set()
            for e in star_arcs:
                arcverts.add(e.source)
                arcverts.add(e.target)
            tmp = neighs.symmetric_difference(arcverts)
            assert len(tmp) == 1
            starvert = tmp.pop()
            star_verts.append(starvert)

    S.vs["star_vert"] = 0
    S.vs[star_verts]["star_vert"] = 1

    S.es.select(_from_in=S.vs.select(star_vert=1)).delete()
    S.es.select(_to_in=S.vs.select(star_vert=1)).delete()

    return S

def main(string_gml, overlap_gml, output_string_gml):
    S = read_graph_gml(string_gml)
    R = read_graph_gml(overlap_gml)

    #  es = len(S.es)

    #  while True:
        #  star_resolution(S, R)
        #  if len(S.es) == es:
            #  break
        #  es = len(S.es)

    S = star_resolution(S, R)
    S.write_gml(output_string_gml)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <string.gml> <overlap.gml> <stars.gml>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    string_gml = sys.argv[1]
    overlap_gml = sys.argv[2]
    output_string_gml = sys.argv[3]

    main(string_gml, overlap_gml, output_string_gml)

    sys.exit(0)
