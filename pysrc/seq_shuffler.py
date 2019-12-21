from Bio import SeqIO
import random
from pathlib import Path

def main():
    # f = '/Users/esaliya/sali/data/cog/uniqs' \
    #     '/1769181_unique_of_1785722_prot2003-2014.fa'

    dir = '/Users/esaliya/sali/data/isolates'
    file = 'isolates_proteins_archaea.fasta'

    f = f"{dir}/{file}"
    fpath = Path(f)
    name = fpath.name
    dir = fpath.parent

    limit = 'inf'
    seqs = list()
    with open(f, "r") as infh:
        seq_count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if limit != 'inf' and seq_count == limit:
                break
            seqs.append(record)
            seq_count += 1

    of = f'{dir}/shuffled_{name}'

    random.shuffle(seqs)
    with open(of, "w") as ofh:
        for record in seqs:
            SeqIO.write(record, ofh, "fasta")


if __name__ == '__main__':
    main()
