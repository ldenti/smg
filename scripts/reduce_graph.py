import sys
from intervaltree import Interval, IntervalTree


def main():
    # assuming a single ref path (e.g., 1 chromosome)
    gfa_path = sys.argv[1]

    k = 0  # FIXME: if we use a value>0, we must be sure to extend to a node in the ref path

    nodes_l = {}
    ref_path_name = ""
    ref_path = []
    tree = IntervalTree()
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            nodes_l[int(idx)] = len(seq)
        elif line.startswith("P"):
            _, tidx, nodes, _ = line.strip("\n").split("\t")
            if not tidx.endswith("_R1"):
                ref_path_name = tidx
                ref_path = [int(x[:-1]) for x in nodes.split(",")]
                continue
            plus = nodes[-1] == "+"
            nodes = [int(x[:-1]) for x in nodes.split(",")]
            if not plus:
                nodes.reverse()
            assert all(b >= a for a, b in zip(nodes[:-1], nodes[1:]))
            tree[nodes[0] - k : nodes[-1] + k + 1] = 1

    # Get reference position for each node on the reference path
    ref_positions = {}
    p = 0
    for n in ref_path:
        ref_positions[n] = p
        p += nodes_l[n]

    print("#genes:", len(tree), file=sys.stderr)
    tree.merge_overlaps()
    i_ = 1
    clusters = {}
    for i in tree:
        clusters[(i.begin, i.end)] = i_
        i_ += 1
    assert len(tree) == len(clusters)
    print("#genic loci:", len(clusters), file=sys.stderr)

    gfa = sys.stdout
    curr_node_idx = -1
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            idx = int(idx)
            curr_node_idx = max(curr_node_idx, idx)
            hit = tree[idx]
            assert len(hit) <= 1
            if len(hit) > 0:
                hit = list(hit)[0]
                c = clusters[(hit.begin, hit.end)]
                line = line[:-1] + f"\tGL:i:{c}"
                if idx in ref_positions:
                    line += f"\tRP:i:{ref_positions[idx]}"
                line += "\n"
                print(line, file=gfa, end="")
        elif line.startswith("L"):
            _, idx1, strand1, idx2, strand2, cigar = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            hit = len(tree[idx1]) > 0 and len(tree[idx2]) > 0
            if hit:
                print(line, file=gfa, end="")
        elif line.startswith("P"):
            _, tidx, nodes, _ = line.strip("\n").split("\t")
            if tidx.endswith("_R1"):
                print(line, file=gfa, end="")


if __name__ == "__main__":
    main()
