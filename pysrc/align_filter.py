import csv

from Bio import SeqIO
import os
from pathlib import Path


def filter_alignments(dir, sub_dir, file, prefix, pid_cut, len_coverage_cut):
    f = f"{dir}/{sub_dir}/{file}"
    of = f"{dir}/{sub_dir}/{prefix}_ani30_lencov70_of_{file}"

    with open(f, "r") as infh, open(of, 'w') as outfh:
        n = 0
        csv_reader = csv.reader(infh, delimiter=',')
        next(csv_reader, None)  # Skip header
        csv_writer = csv.writer(outfh, delimiter=',')
        line_count = 0
        for tup in csv_reader:
            g_col_idx, g_row_idx, pid, col_seq_len, row_seq_len, \
            col_seq_align_len, row_seq_align_len, num_gap_opens, \
            col_seq_len_coverage, row_seq_len_coverage = tup

            pid = float(pid)
            if pid >= pid_cut \
                    and float(row_seq_len_coverage) >= \
                    len_coverage_cut \
                    and float(col_seq_len_coverage) > len_coverage_cut:
                csv_writer.writerow([g_col_idx, g_row_idx, pid/100.0])


def main():
    # The alignment should have an ANI equal to or greater than 30%
    pid_cut = 30
    # The alignment should cover 70% of the length of
    # the both sequences
    len_coverage_cut = 0.7

    dir = '/Users/esaliya/sali/git/github/esaliya/cpp/lbl.distal/pysrc/data' \
          '/cori/scope/distal/w_sub/duprem'
    subs=1600
    alignm='xa'
    sub_dir = f'{alignm}_subs{subs}'
    file = f'{alignm}_shuff_subs{subs}_align.txt'

    filter_alignments(dir, sub_dir, file, "scope", pid_cut,
                      len_coverage_cut)


if __name__ == '__main__':
    main()
