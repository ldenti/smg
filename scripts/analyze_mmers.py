import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    fpath = sys.argv[1]
    mmers = {}
    for line in open(fpath):
        mmer, ccomp = line.strip("\n").split(" ")
        mmer = int(mmer)
        ccomp = int(ccomp)

        if mmer not in mmers:
            mmers[mmer] = [0,set()]
        mmers[mmer][0] += 1
        mmers[mmer][1].add(ccomp)

    data = [] # occurences, number of genic loci
    for v in mmers.values():
        data.append([v[0], len(v[1])])

    data = pd.DataFrame(data, columns = ["N", "G"])
    print(data.describe())
    plt.boxplot(data, vert=True, patch_artist=True, labels=["N", "G"])
    plt.savefig(fpath + ".boxplot.png")
        

if __name__ == "__main__":
    main()
