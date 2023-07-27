
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

def get_linear_chains(G, maxlen):
    chains = []
    for b in G.vs.select(_indegree_ge=3):
        for h in b.successors():
            branch = b.index
            head = h.index
            last = branch
            chain = [branch]
            valid_chain = True
            while h.indegree() == 2:
                if len(chain) >= maxlen-1:
                    valid_chain = False
                    break
                chain.append(h.index)
                u, v = h.successors()
                h = u if u.index != last else v
                last = chain[-1]
            chain.append(h.index)
            if valid_chain:
                chains.append(chain)
    return chains

def alternate_path(G, s, t, mask, maxlen):
    visited = set(mask)
    visited.add(s)
    frontier = [s]
    next_frontier = []
    level = 1
    while len(frontier) > 0:
        if level >= maxlen:
            break
        for u in frontier:
            for vtx in G.vs[u].successors():
                v = vtx.index
                if v == t:
                    return True
                if not v in visited:
                    next_frontier.append(v)
                    visited.add(v)
        frontier = next_frontier[:]
        next_frontier = []
        level += 1
    return False

G = read_graph_gml("G.gml")
G.vs["bubble"] = 0
bubble_pairs = dict()
for chain in get_linear_chains(G, 20):
    if len(chain) <= 2: continue
    u = chain[0]
    v = chain[-1]
    if not u in bubble_pairs:
        bubble_pairs[u] = set()
    if not v in bubble_pairs:
        bubble_pairs[v] = set()
    if v in bubble_pairs[u] and u in bubble_pairs[v]:
        continue
    else:
        bubble_pairs[u].add(v)
        bubble_pairs[v].add(u)
    visited = set(chain[1:-1])
    if alternate_path(G, u, v, visited, 20):
        G.vs[chain]["bubble"] = 1
        print(chain)

G.es.select(_from_in=G.vs.select(bubble=1)).delete()
G.es.select(_to_in=G.vs.select(bubble=1)).delete()

G.write_gml("Z.gml")
