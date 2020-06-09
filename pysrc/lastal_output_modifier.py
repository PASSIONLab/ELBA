from Bio import SeqIO
from pathlib import Path
import csv
import numpy as np
import functools
import operator
import time


def main():
    overlaps_fname = '/Users/esaliya/sali/git/github/esaliya/cpp/lbl.dibella/pysrc/data/cori/scope/mmseqs2.shuff/filtered.gt30pid.gt70lencov.scope.77k.shuff.mmseqs2.out.txt'

    p = Path(overlaps_fname)

    fixed_overlaps_fname = str(p.parent / Path("fixed_" + p.name))

    seqs_fname = '/Users/esaliya/sali/data/scope/uniqs' \
                 '/all/shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa'

    seq_id_2_idx = {}

    with open(seqs_fname, "r") as seqf:
        count = 0
        for record in SeqIO.parse(seqf, "fasta"):
            seq_id_2_idx[record.id] = count
            count += 1

    with open(overlaps_fname, 'rt') as csv_file, \
            open(fixed_overlaps_fname, 'w', newline='') as fixed_csv:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        fixed_csv_writer = csv.writer(fixed_csv, delimiter=' ')
        fixed_csv_writer.writerow(['g_col_idx', 'g_row_idx', 'pid'])
        for tup in csv_reader:
            tup[0] = seq_id_2_idx[tup[0]]
            tup[1] = seq_id_2_idx[tup[1]]
            # fixed_csv_writer.writerow([tup[0], tup[1], tup[2]])
            fixed_csv_writer.writerow(tup)


if __name__ == '__main__':
    main()