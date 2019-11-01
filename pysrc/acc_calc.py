from Bio import SeqIO
import csv
import numpy as np
import functools
import operator
import time


def main():
    overlaps_fname = '/Users/esaliya/sali/git/github/esaliya/cpp/lbl.pisa/pysrc/data/cori/scope/with_sub/k6/subs100/ba2_shuff/ba_shuffled_subs100_align.txt'
    seqs_fname = '/Users/esaliya/sali/data/scope/uniqs' \
                 '/all/77040_unique_of_243813_astral-scopedom-seqres' \
                 '-gd-all-2.07-stable.fa'
    pid_cut = 30
    log_freq = 100000
    align = True

    # All super-families dictionary. Each super family entry will have list,
    # where the first element is the number of families in it and the second
    # is a dictionary for each of its families
    # (family name -> sequence count in family).
    all_sfs = {}
    fam_names = []
    sf_names = []
    limit = 'all'
    count = read_seqs(all_sfs, fam_names, limit, seqs_fname, sf_names)
    tot_pairs_ex_diag = (count * count) - count

    num_sf = len(all_sfs.keys())
    num_fam = sum([len(all_sfs[k][1].keys()) for k in all_sfs])

    print("Total pairs ex diag: ", tot_pairs_ex_diag)
    print("Num super families:  ", num_sf)
    print("Num families:        ", num_fam)

    all_sf_fams_seq_counts = [list(all_sfs[k][1].values()) for k in all_sfs]
    fam_seq_counts = np.array(
        functools.reduce(operator.iconcat, all_sf_fams_seq_counts, []))
    sf_seq_counts = np.array([sum(sf_fams_seq_counts) for sf_fams_seq_counts in
                              all_sf_fams_seq_counts])
    # print(fam_seq_counts)
    # print(sf_seq_counts)

    # Number of family pairs (top triangle only, excludes diagonal as well)
    num_fam_pairs = np.sum(fam_seq_counts * (fam_seq_counts - 1) / 2)
    # Number of super-family pairs
    # (top triangle only, excludes diagonal as well).
    # Includes family pairs too.
    num_sf_pairs = np.sum(sf_seq_counts * (sf_seq_counts - 1) / 2)
    num_sf_only_pairs = num_sf_pairs - num_fam_pairs

    upper_half_pairs_ex_diag = count * (count - 1) / 2

    num_A, num_B, num_C = read_csv_file(fam_names, log_freq, overlaps_fname,
                                        pid_cut, sf_names,
                                        upper_half_pairs_ex_diag, align)
    print()
    print("Output sets ...")
    print("  Set A: ", num_A)
    print("  Set B: ", num_B)
    print("  Set C: ", num_C)

    print()
    print("Upper (U) half ex diag ...")
    print("  All pairs:                             ", upper_half_pairs_ex_diag)
    print("  Family only pairs:                     ", num_fam_pairs)
    print("  Super family (including family) pairs: ", num_sf_pairs)
    print("  Super family only pairs:               ", num_sf_only_pairs)
    print()
    print("  Family only ratio:                     ", num_fam_pairs /
          upper_half_pairs_ex_diag)
    print("  Super family (including family) ratio  ", num_sf_pairs /
          upper_half_pairs_ex_diag)

    # Overall accuracy
    recall = num_A / num_fam_pairs
    precision = (num_A + num_B) / (num_A + num_B + num_C)
    print()
    print("Overall recall:                                   ",
          round(100.0 * recall, 2), "%")
    print("Overall precision:                                ", round(
        100.0 * precision, 2), "%")

    # Family only accuracy
    fam_recall = recall
    fam_precision = num_A / (num_A + num_B + num_C)
    print()
    print("Family only recall:                               ", round(
        100.0 * fam_recall, 2), "%")
    print("Family only precision:                            ", round(
        100.0 * fam_precision, 2), "%")

    # Super family (including family) accuracy
    sf_recall = (num_A + num_B) / (num_fam_pairs + num_sf_only_pairs)
    sf_precision = precision
    print()
    print("Super family (including family) recall:           ", round(
        100.0 * sf_recall, 2), "%")
    print("Super family (including family) recall precision: ",
          round(100.0 * sf_precision, 2), "%")


def read_csv_file(fam_names, log_freq, overlaps_fname, pid_cut, sf_names,
                  upper_half_pairs_ex_diag, align):
    print()
    print("Reading CSV ...")
    t = time.process_time()
    num_A, num_B, num_C = 0, 0, 0
    with open(overlaps_fname, 'rt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        print("  CSV Header\n    ", next(csv_reader))  # ignore header
        line_count = 0
        for tup in csv_reader:

            if align:
                g_col, g_row, pid, col_seq_len, row_seq_len, \
                    col_seq_align_len, row_seq_align_len = tup
            else:
                g_col, g_row, common_kmer_count = tup

            if align and float(pid) < pid_cut:
                continue

            g_col = int(g_col)
            g_row = int(g_row)
            if sf_names[g_col] == sf_names[g_row]:
                if fam_names[g_col] == fam_names[g_row]:
                    num_A += 1
                else:
                    num_B += 1
            else:
                num_C += 1
            line_count += 1
            if line_count % log_freq == 0:
                elapsed = time.process_time() - t
                print("     Lines ", line_count, " of ",
                      upper_half_pairs_ex_diag,
                      " (",
                      round((line_count * 100.0 / upper_half_pairs_ex_diag), 2),
                      "%) took ", round(elapsed, 4), "s")
    print("  Total line count:    ", line_count)
    print("  Total CSV read time: ", round((time.process_time() - t), 2), "s")
    return num_A, num_B, num_C


def read_seqs(all_sfs, fam_names, limit, seqs_fname, sf_names):
    with open(seqs_fname, "r") as seqf:
        count = 0
        for record in SeqIO.parse(seqf, "fasta"):
            if count != 'all' and count == limit:
                break
            l_idx = record.description.index(" ")
            r_idx = record.description.index(" ", l_idx + 1)
            cls, fold, sf, fam = record.description[l_idx: r_idx].split('.')
            fam_names.append(fam)
            sf_names.append(sf)
            if sf in all_sfs:
                sf_fams = all_sfs[sf][1]
                if fam in sf_fams:
                    sf_fams[fam] += 1
                else:
                    sf_fams[fam] = 1
                all_sfs[sf][0] += 1
            else:
                all_sfs[sf] = [1, {fam: 1}]

            count += 1
    print("Read ", count, " sequences")
    return count


if __name__ == '__main__':
    main()
