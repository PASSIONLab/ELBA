#!/usr/bin/env python

import sys
import numpy as np

base_for = "ACGT"
base_rev = "TGCA"
comp_tab = str.maketrans(base_for, base_rev)

def create_data(genome_length, sequencing_depth, avg_read_length, sd_read_length):

    num_reads = int((genome_length * sequencing_depth) / avg_read_length)

    rlengths = np.random.normal(avg_read_length, sd_read_length, num_reads).astype(int)
    rpositions = np.random.randint(0, genome_length - avg_read_length, num_reads)

    rtuples = [] # (read id, genome position, read length)

    for i in range(num_reads):
        rlen = rlengths[i]
        rpos = rpositions[i]
        if rpos + rlen > genome_length:
            rlen = genome_length - rpos
        rtuples.append((i, rpos, rlen))

    genome = "".join("ACGT"[c] for c in np.random.randint(0, 4, genome_length))

    return genome, rtuples

def get_overlaps(read_tuples):

    num_reads = len(read_tuples)
    ordered_read_tuples = sorted(read_tuples, key=lambda x: x[1])

    overlaps = []

    for Q in range(num_reads):
        idQ, posQ, lenQ = ordered_read_tuples[Q]
        for T in range(Q+1, num_reads):
            idT, posT, lenT = ordered_read_tuples[T]
            if posQ + lenQ <= posT:
                break

            begQ = posT - posQ
            begT = 0
            endQ = min(lenQ, posT + lenT - posQ)
            endT = min(posQ + lenQ - posT, lenT)

            overlaps.append((idQ, lenQ, begQ, endQ, idT, lenT, begT, endT))
            overlaps.append((idT, lenT, begT, endT, idQ, lenQ, begQ, endQ))

    return overlaps

def write_overlaps(overlaps, num_reads):
    with open("overlaps.txt", "w") as f:
        for overlap in sorted(overlaps, key=lambda x: num_reads*x[0] + x[3]):
            idQ, lenQ, begQ, endQ, idT, lenT, begT, endT = overlap
            f.write("{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\t{}\n".format(idQ+1, lenQ, begQ, endQ, idT+1, lenT, begT, endT))

def write_data(genome, read_tuples):
    with open("ref.fa", "w") as f:
        f.write(">ref\n{}\n".format(genome))
    with open("reads.fa", "w") as f:
        for rid, rpos, rlen in read_tuples:
            s = genome[rpos:rpos+rlen]
            f.write(">{}\tpos={}\tlen={}\n{}\n".format(rid+1, rpos, rlen, s))

def main():

    genome_length = 100000
    sequencing_depth = 26
    avg_read_length = 8000
    sd_read_length = 1000

    genome, read_tuples = create_data(genome_length, sequencing_depth, avg_read_length, sd_read_length)
    overlaps = get_overlaps(read_tuples)
    num_reads = len(read_tuples)
    write_overlaps(overlaps, num_reads)
    write_data(genome, read_tuples)

main()
