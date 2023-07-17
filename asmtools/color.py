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

S = read_graph_gml("S.gml")

S.es["bad"] = 0

def remove_edge(S, u, v):
    es = S.get_eids([(u, v), (v, u)])
    S.es[es]["bad"] = 1

remove_edge(S, 1720, 1471)
remove_edge(S, 1471, 966)

remove_edge(S, 1063, 8)
remove_edge(S, 8, 1511)

remove_edge(S, 115, 1556)
remove_edge(S, 1556, 1835)

remove_edge(S, 118, 1227)
remove_edge(S, 1227, 115)

remove_edge(S, 1175, 1316)
remove_edge(S, 860, 0)
remove_edge(S, 0, 1175)
remove_edge(S, 96, 1175)

remove_edge(S, 993, 287)
remove_edge(S, 287, 340)
remove_edge(S, 340, 1871)

remove_edge(S, 833, 692)
remove_edge(S, 692, 1073)
remove_edge(S, 1073, 737)


S.write_gml("G.gml")

