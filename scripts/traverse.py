import sys


def main():
    tree_path = sys.argv[1]
    gfa_mmer_path = sys.argv[2]
    fx_mmer_path = sys.argv[3]

    print("Loading tree...", file=sys.stderr)
    comp_to_leaf = {}
    for line in open(tree_path):
        if line.startswith("-"):
            _, leaf, _, components, _ = line.strip("\n").split(" ")
            for c in components.split(","):
                comp_to_leaf[int(c)] = int(leaf)

    print("Loading GFA minimizers...", file=sys.stderr)
    minimizers = {}
    for line in open(gfa_mmer_path):
        mmer, locus = line.strip("\n").split(" ")
        mmer, locus = int(mmer), int(locus)
        leaf = comp_to_leaf[locus]
        minimizers[mmer] = (
            minimizers[mmer] | set([leaf]) if mmer in minimizers else set([leaf])
        )

    for line in open(fx_mmer_path):
        qname, *mmers = line.strip("\n").split(" ")
        if len(mmers) == 0:
            continue
        loci = [
            minimizers[int(mmer)] if int(mmer) in minimizers else set()
            for mmer in mmers
        ]
        print(sum([len(x) for x in loci]))


if __name__ == "__main__":
    main()
