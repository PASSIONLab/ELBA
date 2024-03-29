#!/usr/bin/env python

import sys
import numpy as np

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def create_data(genome_length, sequencing_depth, avg_read_length, sd_read_length):

    num_reads = int((genome_length * sequencing_depth) / avg_read_length)

    rstrands = np.random.randint(0, 2, num_reads)
    rlengths = np.random.normal(avg_read_length, sd_read_length, num_reads).astype(int)
    rpositions = np.random.randint(0, genome_length - avg_read_length, num_reads)

    rtuples = [] # (read id, genome position, read length, read strand)

    for i in range(num_reads):
        rlen = rlengths[i]
        rpos = rpositions[i]
        rstrand = rstrands[i]
        if rpos + rlen > genome_length:
            rlen = genome_length - rpos
        rtuples.append((i, rpos, rlen, rstrand))

    genome = "".join("ACGT"[c] for c in np.random.randint(0, 4, genome_length))

    return genome, rtuples

def get_overlaps(read_tuples):

    num_reads = len(read_tuples)
    ordered_read_tuples = sorted(read_tuples, key=lambda x: x[1])

    overlaps = []

    for Q in range(num_reads):
        idQ, posQ, lenQ, strandQ = ordered_read_tuples[Q]
        for T in range(Q+1, num_reads):
            idT, posT, lenT, strandT = ordered_read_tuples[T]
            if posQ + lenQ <= posT:
                break

            begQ = posT - posQ
            begT = 0
            endQ = min(lenQ, posT + lenT - posQ)
            endT = min(posQ + lenQ - posT, lenT)

            rc = int(strandQ != strandT)

            if begQ <= begT and lenQ - endQ <= lenT - endT:
                direction = -1
                directionT = -1
                suffix = 0
                suffixT = 0
            elif begQ >= begT and lenQ - endQ >= lenT - endT:
                direction = -1
                directionT = -1
                suffix = 0
                suffixT = 0
            elif begQ > begT:
                direction = 0 if rc else 1
                directionT = 0 if rc else 2
                suffix = (lenT - endT) - (lenQ - endQ)
                suffixT = begQ - begT
            else:
                direction = 3 if rc else 2
                directionT = 3 if rc else 1
                suffix = begT - begQ
                suffixT = (lenQ - endQ) - (lenT - endT)

            if strandQ == 1:
                tmp = begQ
                begQ = lenQ - endQ
                endQ = lenQ - tmp

            if strandT == 1:
                tmp = begT
                begT = lenT - endT
                endT = lenT - tmp

            if idQ < idT:
                overlaps.append((idQ, lenQ, begQ, endQ, rc, idT, lenT, begT, endT, direction, suffix))
            else:
                overlaps.append((idT, lenT, begT, endT, rc, idQ, lenQ, begQ, endQ, directionT, suffixT))

    return overlaps

def write_overlaps(overlaps, num_reads):
    with open("overlaps.txt", "w") as f:
        for overlap in sorted(overlaps, key=lambda x: (x[0], x[5])):
            idQ, lenQ, begQ, endQ, rc, idT, lenT, begT, endT, d, s = overlap
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(idQ+1, idT+1, lenQ, begQ, endQ, "+-"[int(rc)], lenT, begT, endT, d, s))

def write_data(genome, read_tuples):
    with open("ref.fa", "w") as f:
        f.write(">ref\n{}\n".format(genome))
    with open("reads.fa", "w") as f:
        for rid, rpos, rlen, rstrand in read_tuples:
            s = genome[rpos:rpos+rlen]
            if rstrand == 1: s = s.translate(comp_tab)[::-1]
            f.write(">{}\tpos={}\tlen={}\tstrand={}\n{}\n".format(rid+1, rpos, rlen, rstrand, s))

def find_contained_reads(overlaps):
    contained = set()
    for idQ, lenQ, begQ, endQ, rc, idT, lenT, begT, endT, _, _ in overlaps:
        if begQ == 0 and endQ == lenQ:
            contained.add(idQ)
        elif begT == 0 and endT == lenT:
            contained.add(idT)
    return contained

def main():

    np.random.seed(313)

    genome_length = 100000
    sequencing_depth = 26
    avg_read_length = 8000
    sd_read_length = 1000

    genome, read_tuples = create_data(genome_length, sequencing_depth, avg_read_length, sd_read_length)
    overlaps = get_overlaps(read_tuples)
    num_reads = len(read_tuples)
    write_overlaps(overlaps, num_reads)
    write_data(genome, read_tuples)

    contained = find_contained_reads(overlaps)
    with open("contained.txt", "w") as f:
        for rid in sorted(contained):
            f.write("{}\n".format(rid+1))

main()
