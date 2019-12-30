from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def main():
    f = '/Users/esaliya/sali/data/cog/prot2003-2014.fa'
    fpath = Path(f)
    dir = fpath.parent
    name_only = fpath.stem
    name = fpath.name

    limit = 'inf'
    seqs = defaultdict(list)
    dups = {}
    with open(f, "r") as infh:
        seq_count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if limit != 'inf' and seq_count == limit:
                break
            seq_upper = str(record.seq.upper())
            seqs[seq_upper].append(record.id)
            if len(seqs[seq_upper]) == 2:
                # Mark the seqs as a duplicate the moment you find
                # out there is another sequence id corresponding to this
                # sequence string, hence the number 2 in the condition.
                # It also avoids adding the sequence once finding it's a
                # duplicate.
                dups[seq_upper] = 1  # 1 is a flag that we'll use later
            seq_count += 1

    # Here, dup_count is just the number of sequences that have multiple
    # copies but not the number of copies itself.
    dup_count = len(dups)

    # This, count, will hold the total number of copies these multi-copy
    # sequences have.
    count = 0
    for dup in dups:
        tmp = len(seqs[dup])
        count += tmp if tmp > 1 else 0

    # Unique sequence count.
    santized_seq_count = (seq_count - count) + dup_count
    dup_count = count - dup_count  # You can keep one copy of each duplicate

    of = f'{dir}/{santized_seq_count}_unique_of_{seq_count}_{name}'
    dup_of = f'{dir}/{dup_count}_duplicates_of_{seq_count}_{name_only}.txt'

    with open(f, "r") as infh, open(of, "w") as ofh, open(dup_of,
                                                          "w") as dup_ofh:
        seq_count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if limit != 'inf' and seq_count == limit:
                break
            seq_upper = str(record.seq.upper())
            if seq_upper not in dups:
                SeqIO.write(record, ofh, "fasta")
            else:
                if dups[seq_upper] == 1:
                    SeqIO.write(record, ofh, "fasta")
                    dups[seq_upper] = 0
                else:
                    continue
            seq_count += 1

        dup_ofh.write(f"Duplicate Count: {dup_count}\n\n")
        for seq_upper, ids in seqs.items():
            if len(ids) > 1:
                dup_ofh.write(ids[0]+":")
                dup_ofh.write(",".join(ids[1:]))
                dup_ofh.write("\n" + str(seq_upper) + "\n")


if __name__ == '__main__':
    main()
