from igraph import *
import pdb

def asg_to_igraph(asg, seqs, names):
    n = len(asg)
    G = Graph(n, directed=True)
    es, ds, ls = [], [], []
    G.vs["readname"] = names
    G.vs["readlen"] = [len(s) for s in seqs]
    for u in asg:
        for v in asg[u]:
            es.append((u,v))
            ds.append(asg[u][v][0])
            ls.append(asg[u][v][1])
    G.add_edges(es)
    G.es["direction"] = ds
    G.es["length"] = ls
    return G

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def revcomp(s): return s.translate(comp_tab)[::-1]

def read_paf(paf_fname):
    overlaps = []
    types = (str, int, int, int, lambda x: int(x[0] == '-'), str, int, int, int)
    for line in open(paf_fname, "r"):
        tokens = line.rstrip().split()[:9]
        overlaps.append(tuple(t(tok) for t, tok in zip(types, tokens)))
    return overlaps

def read_fasta(fasta_fname):
    seqs, names, seq, name = [], [], [], ""
    for line in open(fasta_fname):
        if line[0] == ">":
            if len(seq) > 0:
                seqs.append("".join(seq))
                names.append(name)
            name = line.lstrip(">").split()[0].rstrip()
            seq = []
        else: seq.append(line.rstrip())
    if len(seq) > 0:
        seqs.append("".join(seq))
        names.append(name)
    return seqs, names

def asg_delete_verts(asg, vdel):
    g = {i : {} for i in range(len(asg))}
    for u in asg:
        if u in vdel: continue
        for v in asg[u]:
            if v in vdel: continue
            g[u][v] = asg[u][v]
    return g

def gen_asm_graph(overlaps, seqs, names):
    n = len(seqs)
    namemap = {names[i] : i for i in range(n)}
    asg = {i : {} for i in range(n)}
    vdel = set()
    for overlap in overlaps:
        qn, ql, qb, qe, rc, tn, tl, tb, te = overlap
        if qn == tn: continue
        qi, ti = namemap[qn], namemap[tn]
        assert ql == len(seqs[qi]) and tl == len(seqs[ti])
        b1, e1, l1 = qb, qe, ql
        b2, e2, l2 = (tb, te, tl) if rc == 0 else (tl - te, tl - tb, tl)
        if b1 <= b2 and l1 - e1 <= l2 - e2:
            vdel.add(qi)
            continue
        elif b1 >= b2 and l1 - e1 >= l2 - e2:
            vdel.add(ti)
            continue
        elif b1 > b2:
            asg[qi][ti] = (1 if rc == 0 else 0, b1 - b2)
            asg[ti][qi] = (2 if rc == 0 else 0, (l2 - e2) - (l1 - e1))
        else:
            asg[ti][qi] = (1 if rc == 0 else 3, b2 - b1)
            asg[qi][ti] = (2 if rc == 0 else 3, (l1 - e1) - (l2 - e2))

    return asg_delete_verts(asg, vdel)

def get_root(asg, x, z, visited):
    ll = 0
    while True:
        visited.add(x)
        ll += 1
        if len(asg[x]) != 2:
            if len(asg[x]) == 1:
                z = x
            break
        z = x
        x = [k for k in asg[x] if not k in visited][0]
    return z, ll

def get_closest_root(asg, v):
    assert len(asg[v]) <= 2
    a = [j for j in asg[v] if len(asg[j]) <= 2]
    if len(a) <= 1: return v
    x, y = a
    visited = set([v])

    lroot, ll = get_root(asg, x, v, visited)
    rroot, lr = get_root(asg, y, v, visited)

    return lroot if ll < lr else rroot

def get_utg_comps(asg):
    n = len(asg)
    mark = [False] * n
    utgs = []
    for v in range(n):
        if mark[v] or len(asg[v]) == 0 or len(asg[v]) >= 3: continue
        utg = []
        r = get_closest_root(asg, v)
        while True:
            mark[r] = True
            utg.append(r)
            if len(asg[r]) != 2 and len(utg) > 1: break
            lr = [x for x in asg[r] if not mark[x] and len(asg[x]) <= 2]
            if len(lr) != 1:
                assert len(lr) == 0
                break
            r = lr[0]
        utgs.append(utg)
    return utgs

def asm_utg(asg, utg, seqs):
    if len(utg) == 1: return seqs[utg[0]], [(utg[0], 0, 0, 0, len(seqs[utg[0]]))]
    p = []
    q = len(utg)
    ld = None
    gfainfo = []
    upos = 0
    for i in range(q-1):
        u, v = utg[i], utg[i+1]
        d, l = asg[u][v]
        rc = int((d>>1)&1)
        s = seqs[u]
        if rc == 1: s = revcomp(s)
        p.append(s[:l])
        ld = d
        gfainfo.append((u, upos, rc, 0, l))
        upos += l
    s = seqs[utg[-1]]
    rc = 1 - int(ld&1)
    if rc == 1: s = revcomp(s)
    p.append(s)
    gfainfo.append((utg[-1], upos, rc, 0, len(s)))
    return "".join(p), gfainfo

def gen_ug(asg, seqs, names):
    ug = dict()
    bmap = dict() # maps graph indices of one-read utgs to their unitig graph indices
    umap = dict() # maps unitig graph indices to linear list of read indices
    utgs = get_utg_comps(asg)
    utgseqs = []
    gfainfos = []
    i = 0
    for utg in utgs:
        l, r = utg[0], utg[-1]
        if l == r: bs = [x for x in asg[l]] + [l]
        else: bs = [x for x in asg[l] if x != utg[1]] + [x for x in asg[r] if x != utg[-2]]
        for b in bs:
            if not b in bmap:
                bmap[b] = i
                ug[i] = dict()
                us, gfainfo = asm_utg(asg, utg, seqs)
                utgseqs.append(us)
                gfainfos.append(gfainfo)
                i += 1
        if l == r:
            for d in asg[l]:
                ug[bmap[l]][bmap[d]] = asg[l][d]
                ug[bmap[d]][bmap[l]] = asg[d][l]
        else:
            umap[i] = utg
            ug[i] = dict()
            us, gfainfo = asm_utg(asg, utg, seqs)
            utgseqs.append(us)
            gfainfos.append(gfainfo)
            i += 1
            for d in asg[l]:
                if d == utg[1]: continue
                ug[bmap[d]][i-1] = (asg[d][l][0] | 1, len(seqs[d]) - asg[d][l][1])
                ug[i-1][bmap[d]] = (asg[l][d][0] | 2, len(seqs[d]) - asg[d][l][1])
            for d in asg[r]:
                if d == utg[-2]: continue
                ug[i-1][bmap[d]] = (asg[r][d][0] & 1, len(seqs[d]) - asg[d][r][1])
                ug[bmap[d]][i-1] = (asg[d][r][0] & 2, len(seqs[d]) - asg[d][r][1])
    return ug, bmap, umap, utgseqs, gfainfos

def ug2gfa(gfa_fname, ug, bmap, umap, utgseqs, gfainfos, names):
    with open(gfa_fname, "w") as f:
        for i, s in enumerate(utgseqs):
            f.write("S\tutg{}\t{}\tLN:i:{}\n".format(i+1, s, len(s)))
            for readid, upos, rc, start, end in gfainfos[i]:
                f.write("A\tutg{}\t{}\t{}\t{}\t{}\t{}\n".format(i+1, upos, "+-"[rc], names[readid], start, end))
        for u in ug:
            for v in ug[u]:
                d, l = ug[u][v]
                f.write("L\tutg{}\t{}\tutg{}\t{}\t{}M\n".format(u+1, "+-"[(d>>1)&1], v+1, "-+"[d&1], l))


fasta_fname = "../reads.fa"
paf_fname = "elba.string.paf"
seqs, names = read_fasta(fasta_fname)
overlaps = read_paf(paf_fname)
asg = gen_asm_graph(overlaps, seqs, names)
g = asg_to_igraph(asg, seqs, names)
ug, bmap, umap, utgseqs, gfainfos = gen_ug(asg, seqs, names)

g.vs["comp"] = 0
for k in bmap:
    g.vs[k]["comp"] = bmap[k]+1
for k in umap:
    g.vs[umap[k]]["comp"] = k+1

g.write_gml("g.gml")

ug2gfa("test.gfa", ug, bmap, umap, utgseqs, gfainfos, names)




