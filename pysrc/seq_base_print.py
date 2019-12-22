from Bio import SeqIO
import os
from pathlib import Path


def print_bases(dir, fname, start, count):
    alph = {'A','R','N','D','C','Q','E','G','H','I','L','K','M',
                'F','P','S','T','W','Y','V','B','Z','X','*','J'}
    f = f"{dir}/{fname}"

    with open(f, "r") as infh:
        n = 0
        limit = start + count
        gen = SeqIO.parse(infh, "fasta")
        bases = set()
        # rottens = set()
        for record in gen:
            if n < start:
                n += 1
                continue
            if n < limit:
                for c in record:
                    bases.add(c)
                    if c not in alph:
                        print(n)
                n += 1
            else:
                break

    # print(sorted(rottens))
    print(len(bases), ":", bases)


def main():
    # dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    # file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    # dir = '/Users/esaliya/sali/data/scope/uniqs/all'
    # file = 'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2' \
    #        '.07-stable.fa'

    dir = '/Users/esaliya/sali/data/isolates/archaea'
    file = '2729008_len_lte_2000_in_shuffled_isolates_proteins_archaea.fasta'
    start = 0
    count = 2729008
    print_bases(dir, file, start, count)

    # alph = {'A','R','N','D','C','Q','E','G','H','I','L','K','M',
    #             'F','P','S','T','W','Y','V','B','Z','X','*','J'};

    # bases = {'X', 'I', 'T', 'N', 'C', 'O', 'M', 'F', 'B', 'K', 'D', 'A',
    #          'H', 'Q', 'P', 'L', 'G', 'Y', 'S', 'W', 'Z', 'J', 'R', 'U', 'V', 'E'}

    # print (alph - bases)
    # print (bases - alph)


if __name__ == '__main__':
    main()
