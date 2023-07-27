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

def bubble_walk(G, branch, head):
    b = G.vs[branch]
    h = G.vs[head]
    assert b.indegree() >= 3
    last = branch
    chain = [branch]
    while h.indegree() == 2:
        chain.append(h.index)
        u, v = h.successors()
        h = u if u.index != last else v
        last = chain[-1]
    chain.append(h.index)
    return chain

def label_short_linear_walks(G, maxlen):
    G.vs["short_walk"] = 0
    for b in G.vs.select(_indegree_ge=3):
        for h in b.successors():
            walk = bubble_walk(G, b.index, h.index)
            if len(walk) >= 3 and len(walk) <= maxlen:
                G.vs[walk[1:-1]]["short_walk"] = len(walk)
    return G

def bubble_walks(G, branch):
    b = G.vs[branch]
    assert b.indegree() >= 3
    walks = {}
    for h in b.successors():
        walks[h.index] = bubble_walk(G, branch, h.index)
    return walks

def identify_bubbles(G, maxlen):
    bubbles = []
    visited = set()
    for b in G.vs.select(_indegree_eq=3):
        branch = b.index
        walks = bubble_walks(G, branch)
        info = []
        for head, walk in walks.items():
            if len(walk) >= 3 and len(walk) <= maxlen:
                info.append((head, walk[-1]))
        for i in range(len(info)):
            for j in range(i):
                if info[i][1] == info[j][1]:
                    r1 = walks[info[i][0]]
                    r2 = walks[info[j][0]]
                    if not r1[0] in visited and not r1[-1] in visited:
                        visited.add(r1[0])
                        visited.add(r1[-1])
                        bubbles.append((r1, r2))
    return bubbles

def remove_edge(S, u, v):
    es = S.get_eids([(u, v), (v, u)])
    S.es[es].delete()

def main(input_gml, output_gml, maxlen):

    G = read_graph_gml(input_gml)
    bubbles = identify_bubbles(G, maxlen)

    G.vs["bubble"] = 0
    for bubble in bubbles:
        G.vs[bubble[0]]["bubble"] = 1
        G.vs[bubble[1]]["bubble"] = 2

    G.write_gml(output_gml)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("usage: {} <graph.gml> <graph.nobubbles.gml> <maxlen>\n".format(sys.argv[0]))
        sys.stderr.flush()
        sys.exit(-1)

    input_gml = sys.argv[1]
    output_gml = sys.argv[2]
    maxlen = int(sys.argv[3])

    main(input_gml, output_gml, maxlen)

    sys.exit(0)
