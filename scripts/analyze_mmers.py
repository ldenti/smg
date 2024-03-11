import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    fpath = sys.argv[1]
    mmers = {}
    loci = {}
    for line in open(fpath):
        mmer, ccomp = line.strip("\n").split(" ")
        mmer = int(mmer)
        ccomp = int(ccomp)

        if mmer not in mmers:
            mmers[mmer] = [0,set()]
        mmers[mmer][0] += 1
        mmers[mmer][1].add(ccomp)

        if ccomp not in loci:
            loci[ccomp] = set()
        loci[ccomp].add(mmer)

    data = [] # occurences, number of genic loci
    for v in mmers.values():
        data.append([v[0], len(v[1])])

    data = pd.DataFrame(data, columns = ["N", "G"])
    print(data.describe())
    plt.boxplot(data, vert=True, patch_artist=True, labels=["N", "G"])
    plt.savefig(fpath + ".boxplot.png")

    big_matrix = [[0 for _ in range(len(loci))] for _ in range(len(loci))]
    for i, (l1, mm1) in enumerate(loci.items()):
        for j, (l2, mm2) in enumerate(loci.items()):
            if j < i:
                big_matrix[i][j] = 0
            elif j == i:
                big_matrix[i][j] = 1
            else:
                big_matrix[i][j] = len(mm1 & mm2) / len(mm1 | mm2)
    sns.heatmap(big_matrix)
    plt.savefig(fpath + ".heatmap.png")

if __name__ == "__main__":
    main()
