from Bio import SeqIO
import os
from pathlib import Path


def rem_bad_seqs(dir, fname, bad_seqs):
    f = f"{dir}/{fname}"
    bad_f = f"{dir}/{bad_seqs}"

    rottens = set()
    with open(bad_f, "r") as badfh:
        for line in badfh:
            rottens.add(int(line))

    of = f"{dir}/sanitized_{fname}"

    n = 0
    with open(f, "r") as infh, open(of, 'w') as outfh:
        gen = SeqIO.parse(infh, "fasta")
        for record in gen:
            if n not in rottens:
                SeqIO.write(record, outfh, "fasta")
            n += 1

    ofnew = f"{dir}/sanitized_{n - len(rottens)}_{fname}"
    os.rename(of, ofnew)

def main():
    # dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    # file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    # dir = '/Users/esaliya/sali/data/scope/uniqs/all'
    # file = 'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2' \
    #        '.07-stable.fa'

    dir = '/Users/esaliya/sali/data/isolates/archaea'
    file = 'impure_2729008_len_lte_2000_in_shuffled_isolates_proteins_archaea' \
           '.fasta'
    bad_seq_file = 'bad_seqs.txt'

    rem_bad_seqs(dir, file, bad_seq_file)


if __name__ == '__main__':
    main()
