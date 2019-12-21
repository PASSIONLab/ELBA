from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


def main():
    # dir = '/Users/esaliya/sali/data/cog/uniqs/shuffled'
    # file = 'shuffled_1769181_unique_of_1785722_prot2003-2014.fa'

    dir = '/Users/esaliya/sali/data/isolates'
    file = 'isolates_proteins_archaea.fasta'

    f = f"{dir}/{file}"

    limit = 'inf'

    lengths = list()
    with open(f, "r") as infh:
        count = 0
        for record in SeqIO.parse(infh, "fasta"):
            if limit != 'inf' and count == limit:
                break
            lengths.append(len(record))
            count += 1

    print("INFO: Sequence count: ", count)

    tag = str(count)+"_of_" if limit != 'inf' else "all"
    of = f'{dir}/hist_of_{tag}_{Path(file).stem}.jpg'

    lengths = np.array(lengths)
    plot_sequence_lengths_histogram(lengths, 'auto', of)


def plot_sequence_lengths_histogram(lengths, nbins, fname):
    plt.rcParams.update({'font.size': 30})
    plt.figure(figsize=(50, 20))
    n, bins, patches = plt.hist(lengths, bins=nbins, edgecolor='black', log=False)
    patch_width = bins[1] - bins[0]
    plt.title("Histogram of Sequence Lengths")
    plt.ylim(0.1, np.max(n) + 20)
    # plt.minorticks_on()
    plt.xlabel("Sequence Length")
    plt.ylabel("# Sequences")
    # ax = plt.gca()
    # ax.yaxis.grid(which='both')
    # ax.set_xticks(np.arange(bins[0]+(patch_width/2), bins[-1]+0.5, patch_width*1))
    # ax.set_xticklabels(np.stack((np.around(bins[0:len(bins) - 1], 2),
    #                              np.around(bins[0:len(bins) - 1] + patch_width, 2)), 0).T)
    # for i in range(len(n)):
    #     plt.text(bins[i] + 0.2, n[i] + np.log(2), str(int(n[i])))

    plt.savefig(fname)


if __name__ == '__main__':
    main()

