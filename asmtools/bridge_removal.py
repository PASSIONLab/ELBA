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

def branch_walk(G, branch, head):
    b = G.vs[branch]
    h = G.vs[head]
    assert b.indegree() >= 3
    last = branch
    chain = []
    while h.indegree() == 2:
        chain.append(h.index)
        u, v = h.successors()
        h = u if u.index != last else v
        last = chain[-1]
    return chain

def branch_walks(G, branch):
    b = G.vs[branch]
    assert b.indegree() >= 3
    walks = {}
    for h in b.successors():
        walks[h.index] = len(branch_walk(G, branch, h.index))
    return walks

def get_bridges(G):
    visited = set()
    bridges = set()
    for triple in G.vs.select(_indegree_eq=3):
        for u in triple.successors():
            if not u.index in visited:
                visited.add(u.index)
            else:
                bridges.add(u.index)
    return bridges

def get_lonely_bridges(G, l):
    bridges = get_bridges(G)
    lonely_bridges = []
    for bridge in bridges:
        b = G.vs[bridge]
        if b.indegree() != 2:
            continue
        u, v = b.successors()
        u_walk = branch_walks(G, u.index)
        v_walk = branch_walks(G, v.index)
        ucnt, vcnt = 0, 0
        for _, cnt in u_walk.items():
            if cnt >= l: ucnt += 1
        for _, cnt in v_walk.items():
            if cnt >= l: vcnt += 1
        if ucnt == 2 and vcnt == 2:
            lonely_bridges.append(bridge)
    return lonely_bridges

def remove_edge(S, u, v):
    es = S.get_eids([(u, v), (v, u)])
    S.es[es].delete()

def main(input_gml, output_gml, walklen):

    G = read_graph_gml(input_gml)
    bridges = get_lonely_bridges(G, walklen)

    for bridge in bridges:
        b = G.vs[bridge]
        u,v = [w.index for w in b.successors()]
        remove_edge(G, bridge, u)
        remove_edge(G, bridge, v)

    G.write_gml(output_gml)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <graph.gml> <graph.nobridges.gml> <walklen>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    input_gml = sys.argv[1]
    output_gml = sys.argv[2]
    walklen = int(sys.argv[3])

    main(input_gml, output_gml, walklen)

    sys.exit(0)
