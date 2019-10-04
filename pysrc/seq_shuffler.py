from Bio import SeqIO
import random


def main():
    f = '/Users/esaliya/sali/data/scope/uniqs/all' \
        '/77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa'

    limit = 'inf'
    seqs = list()
    with open(f, "r") as infh:
        seq_count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if limit != 'inf' and seq_count == limit:
                break
            seqs.append(record)
            seq_count += 1

    of = f = '/Users/esaliya/sali/data/scope/uniqs/all/' \
        'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07' \
             '-stable.fa'

    random.shuffle(seqs)
    with open(of, "w") as ofh:
        for record in seqs:
            SeqIO.write(record, ofh, "fasta")


if __name__ == '__main__':
    main()
