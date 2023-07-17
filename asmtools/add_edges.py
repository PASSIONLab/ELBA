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

def get_contig_labelled_graph(S):
    G = S.copy()
    G.vs["contig"] = -1
    G.vs["branch"] = 0
    G.vs["contigroot"] = 0
    H = G.copy()
    branches = H.vs.select(_indegree_ge=3)
    G.vs[branches.indices]["branch"] = 1
    branch_incident = H.es.select(_incident=branches)
    branch_incident.delete()
    G.vs[H.vs.select(_indegree_eq=1).indices]["contigroot"] = 1
    comps = H.components()
    large_contig_membership = Counter(comps.membership).most_common()
    smallid = 0
    for largeid, cnt in large_contig_membership:
        if cnt < 2: break
        G.vs[comps[largeid]]["contig"] = smallid
        smallid += 1
    return G

S = read_graph_gml("S.gml")
S = get_contig_labelled_graph(S)

R = read_graph_gml("R.gml")
edges = R.es.select(_within=G.vs.select(contigroot=1).indices)

S.es["new"] = 0
S.add_edges([e.tuple for e in edges])
S.es.select(new=None)["new"] = 1

S.write_gml("example.gml")

#  def get_contig_labelled_graph(S):
    #  G = S.copy()
    #  G.vs["contig"] = -1
    #  G.vs["branch"] = 0
    #  G.vs["contigroot"] = 0
    #  H = G.copy()
    #  branches = H.vs.select(_indegree_ge=3)
    #  G.vs[branches.indices]["branch"] = 1
    #  branch_incident = H.es.select(_incident=branches)
    #  branch_incident.delete()
    #  G.vs[H.vs.select(_indegree_eq=1).indices]["contigroot"] = 1
    #  comps = H.components()
    #  large_contig_membership = Counter(comps.membership).most_common()
    #  smallid = 0
    #  for largeid, cnt in large_contig_membership:
        #  if cnt < 2: break
        #  G.vs[comps[largeid]]["contig"] = smallid
        #  smallid += 1
    #  return G


