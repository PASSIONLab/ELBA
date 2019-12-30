from Bio import SeqIO
import os
from pathlib import Path


def extract_families(dir, fname):
    f = f"{dir}/{fname}"
    of = f"{dir}/family_info_{Path(fname).stem}.txt"

    with open(f, "r") as seqf, open(of, 'w') as outfh:
        count = 0
        for record in SeqIO.parse(seqf, "fasta"):
            l_idx = record.description.index(" ")
            r_idx = record.description.index(" ", l_idx + 1)
            cls, fold, sf, fam = record.description[l_idx: r_idx].split('.')

            outfh.write(f"{count} {cls}.{fold}.{sf}.{fam}\n");
            count += 1
    print("Read ", count, " sequences")


def main():
    # dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    # file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    dir = '/Users/esaliya/sali/data/scope/uniqs/all'
    file = 'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa'

    # dir = '/Users/esaliya/sali/data/isolates/archaea'
    # file = '2729008_len_lte_2000_in_shuffled_isolates_proteins_archaea.fasta'

    extract_families(dir, file)


if __name__ == '__main__':
    main()
