import sys

def read_fasta(fname):
    seqs = []
    seq = []
    name = ""
    namemap = {}
    namelist = []
    idx = 0

    for line in open(fname, "r"):
        if line[0] == ">":
            if len(seq) > 0:
                seqs.append("".join(seq))
                namemap[name] = idx
                namelist.append(name)
                seq = []
                idx += 1
            name = line.lstrip(">").split()[0].rstrip()
        else: seq.append(line.rstrip())

    if len(seq) > 0:
        seqs.append("".join(seq))
        namemap[name] = idx
        namelist.append(name)

    return seqs, namemap, namelist

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def revcomp(s): return s.translate(comp_tab)[::-1]

correct = 0
incorrect = 0

fname = sys.argv[1]
ksize = int(sys.argv[2])
seqs, namemap, namelist = read_fasta(fname)

with open("B.mtx", "r") as f:
    next(f)
    next(f)
    for line in f.readlines():
        tokens = tuple(int(v) for v in line.rstrip().split())
        idxQ = tokens[0]
        idxT = tokens[1]
        if idxQ == idxT: continue
        seqQ = seqs[idxQ]
        seqT = seqs[idxT]
        begQs = (tokens[2], tokens[4])
        begTs = (tokens[3], tokens[5])
        seedQs = tuple(seqQ[b:b+ksize] for b in begQs)
        seedTs = tuple(seqT[b:b+ksize] for b in begTs)

        if seedQs[0] != seedTs[0] and seedQs[0] != revcomp(seedTs[0]):
            incorrect += 1
            sys.stdout.write("{}\t{}\t{}\t{}\tno_match\n".format(idxQ, idxT, begQs[0], begTs[0]))
        else: correct += 1

        if seedQs[1] != seedTs[1] and seedQs[1] != revcomp(seedTs[1]):
            incorrect += 1
            sys.stdout.write("{}\t{}\t{}\t{}\tno_match\n".format(idxQ, idxT, begQs[1], begTs[1]))
        else: correct += 1

print("correct={}\nincorrect={}\n".format(correct, incorrect))
