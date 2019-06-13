from Bio import SeqIO


def main():
    f = '/Users/esaliya/sali/data/scope/astral-scopedom-seqres-gd-all-2.07' \
        '-stable.fa'
    limit = 1000
    of = f'/Users/esaliya/sali/data/scope/{limit}_astral-scopedom-seqres-gd' \
        f'-all-2.07-stable.fa'

    with open(f, "r") as infh, open(of, 'w') as outfh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if count == limit:
                break
            SeqIO.write(record, outfh, "fasta")
            count += 1


if __name__ == '__main__':
    main()
