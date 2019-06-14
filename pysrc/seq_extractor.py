from Bio import SeqIO


def main():
    dir = '/Users/esaliya/sali/data/scope/uniqs'
    file = '77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa'
    f = f"{dir}/{file}"

    limit = 1000
    of = f"{dir}/{limit}_of_{file}"

    with open(f, "r") as infh, open(of, 'w') as outfh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if count == limit:
                break
            SeqIO.write(record, outfh, "fasta")
            count += 1


if __name__ == '__main__':
    main()
