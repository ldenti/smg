import sys
import copy

import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt


class Node:
    def __init__(self, idx):
        self.idx = idx
        self.childs = [None, None]

    def set_child(self, i, child):
        self.childs[i] = child

    def __repr__(self):

        return f"{self.idx}|{self.childs[0].idx if self.childs[0] != None else -1}|{self.childs[1].idx if self.childs[1] != None else -1}"


def build_tree(T, node):
    if node.childs[0] != None:
        T.add_edge(node.idx, node.childs[0].idx)
        build_tree(T, node.childs[0])
    if node.childs[1] != None:
        T.add_edge(node.idx, node.childs[1].idx)
        build_tree(T, node.childs[1])


def print_tree(node):
    print(node)
    if node.childs[0] != None:
        print_tree(node.childs[0])
    if node.childs[1] != None:
        print_tree(node.childs[1])


def jaccard(S1, S2):
    return len(S1 & S2) / len(S1 | S2)


def main():
    fpath = sys.argv[1]
    gfa_path = sys.argv[2]
    # filter minimizers with more than this occurrences
    w = 5
    # merge clusters if jaccard is greater than this
    s = 0
    # maximum size for clusters (as total number of characters)
    max_cluster_size = 100000

    # Load components sizes from GFA (as total number of characters)
    print("Load GFA...", file=sys.stderr)
    components = {}
    for line in open(gfa_path):
        if line.startswith("S"):
            # FIXME: this is hardcoded and breaks if used on non-pruned GFA
            _, idx, seq, gl, *_ = line.strip("\n").split("\t")
            gl = int(gl.split(":")[2])
            components[gl] = components[gl] + len(seq) if gl in components else len(seq)

    # Load minimizers and loci information
    print("Load minimizers...", file=sys.stderr)
    minimizers = {}
    loci = {}
    for line in open(fpath):
        mmer, locus = line.strip("\n").split(" ")
        mmer = int(mmer)
        locus = int(locus)

        minimizers[mmer] = minimizers[mmer] + 1 if mmer in minimizers else 1
        loci[locus] = loci[locus] | set([mmer]) if locus in loci else set([mmer])

    # Remove overexpressed mmers
    print("Filter minimizers...", file=sys.stderr)
    for locus, mmers in loci.items():
        loci[locus] = mmers - set([mmer for mmer in mmers if minimizers[mmer] > w])

    # Build clusters based on jaccard similarity
    print("Clustering based on jaccard similarity...", file=sys.stderr)
    clusters = []
    used_flags = [0 for _ in loci]
    for i, (l1, mm1) in enumerate(loci.items()):
        if used_flags[l1 - 1]:
            continue
        cluster = [l1]
        used_flags[l1 - 1] = 1
        for j, (l2, mm2) in enumerate(loci.items()):
            if j > i and not used_flags[l2 - 1] and jaccard(mm1, mm2) > 0:
                cluster.append(l2)
                used_flags[l2 - 1] = 1
        clusters.append(cluster)
    clusters.sort(reverse=True, key=lambda x: sum([components[c] for c in x]))

    # Merge clusters based on their size to have the most balanced clusters we can have
    print("Merging clusters...", file=sys.stderr)
    new_clusters = [clusters[0]]
    prev_cluster_size = sum([components[c] for c in clusters[0]])
    for cluster in clusters[1:]:
        # Alternative 1
        # cluster_size = sum([components[c] for c in cluster])
        # if prev_cluster_size + cluster_size > max_cluster_size:
        #     new_clusters.append([])
        #     prev_cluster_size = 0
        # new_clusters[-1].extend(cluster)
        # prev_cluster_size += cluster_size

        # Alternative 2
        if prev_cluster_size > max_cluster_size:
            new_clusters.append([])
            prev_cluster_size = 0
        new_clusters[-1].extend(cluster)
        prev_cluster_size += sum([components[c] for c in cluster])

    # Create new clusters by minimizers merging (useless now)
    # clusters = []
    for i, cluster in enumerate(new_clusters):
        print(
            "-",
            i,
            len(cluster),
            ",".join([str(c) for c in cluster]),
            sum([components[c] for c in cluster]),
        )
        # clusters.append([cluster, set.union(*[loci[c] for c in cluster])])

    # Build tree structure
    print("Building tree...", file=sys.stderr)
    nodes = [Node(str(x)) for x in list(range(len(new_clusters)))]
    new_nodes = []
    while len(nodes) > 1:
        while len(nodes) != 0:
            node1 = nodes.pop(0)
            if len(nodes) == 0:
                new_nodes.append(node1)
            else:
                node2 = nodes.pop(0)
                root = Node(f"{node1.idx},{node2.idx}")
                root.set_child(0, node1)
                root.set_child(1, node2)
                new_nodes.append(root)
        nodes = new_nodes
        new_nodes = []

    # Print tree to stdout
    print("Printing tree...", file=sys.stderr)
    print_tree(root)

    # Draw tree
    print("Drawing tree...", file=sys.stderr)
    root = nodes[0]
    T = nx.DiGraph()
    build_tree(T, root)
    pos = graphviz_layout(T, prog="dot")
    nx.draw(T, pos, with_labels=True)
    plt.savefig(fpath + ".tree.png")


if __name__ == "__main__":
    main()
