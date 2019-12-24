from Bio import SeqIO
import os
from pathlib import Path
import csv


def generate_index(dir, fname):
    f = f"{dir}/{fname}"
    of = f"{dir}/index_of_{Path(fname).stem}.txt"

    with open(f, "r") as infh, open(of, 'w') as outfh:
        csv_writer = csv.writer(outfh, delimiter=' ')
        n = 0
        gen = SeqIO.parse(infh, "fasta")
        for record in gen:
            csv_writer.writerow([n, record.id])
            n += 1


def main():
    # dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    # file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    dir = '/Users/esaliya/sali/data/scope/uniqs/all'
    file = 'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2' \
           '.07-stable.fa'

    # dir = '/Users/esaliya/sali/data/isolates/archaea'
    # file = '2729008_len_lte_2000_in_shuffled_isolates_proteins_archaea.fasta'

    limit = 10000
    cut = 1000

    generate_index(dir, file)


if __name__ == '__main__':
    main()
