from Bio import SeqIO
import os
from pathlib import Path


def extract_by_count(dir, fname, limit):
    f = f"{dir}/{fname}"
    of = f"{dir}/{limit}_of_{fname}"

    with open(f, "r") as infh, open(of, 'w') as outfh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if count == limit:
                break
            SeqIO.write(record, outfh, "fasta")
            count += 1


def extract_by_length(dir, fname, cut):
    f = f"{dir}/{fname}"
    of = f"{dir}/len_lte_{cut}_in_{fname}"

    with open(f, "r") as infh, open(of, 'w') as outfh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if len(str(record.seq)) <= cut:
                SeqIO.write(record, outfh, "fasta")
                count += 1

    ofnew = f"{dir}/{count}_{Path(of).name}"
    outfh.close()
    os.rename(of, ofnew)


def main():
    dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    limit = 10
    cut = 1000

    # extract_by_count(dir, file, limit)
    extract_by_length(dir, file, cut)


if __name__ == '__main__':
    main()
