from Bio import SeqIO


def main():
    dir = '/Users/esaliya/sali/data/scope/uniqs/all'
    file = 'shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa'
    f = f"{dir}/{file}"

    id1 = 288
    id2 = 255

    with open(f, "r") as infh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if count == id1 or count == id2:
                print(count)
                print(record)
                print(record._seq)
            count += 1


if __name__ == '__main__':
    main()