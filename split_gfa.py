import sys


def main():
    gfa_path = sys.argv[1]
    tree_path = sys.argv[2]
    wd = sys.argv[3]

    subgraphs = {}
    component_to_subgraph = {}
    leafs = {}
    subgraph_idx = 0
    for line in open(tree_path):
        if line.startswith("-"):
            _, leaf, size, components, length = line.strip("\n").split(" ")
            leafs[int(leaf)] = [int(x) for x in components.split(",")]
        else:
            node, child1, child2 = line.strip("\n").split("|")
            if child1 != "-1" or child2 != "-1":
                continue
            subgraphs[subgraph_idx] = leafs[int(node)]
            for x in leafs[int(node)]:
                component_to_subgraph[x] = subgraph_idx
            subgraph_idx += 1
    print(subgraphs)

    vertex_to_subgraph = {}
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq, gl, *_ = line.strip("\n").split("\t")
            gl = int(gl.split(":")[2])
            subgraph = component_to_subgraph[gl]
            vertex_to_subgraph[idx] = subgraph

    ofiles = []
    for k in subgraphs:
        ofiles.append(open(f"{wd}/leaf_{k}.gfa", "w"))
    for LINE in open(gfa_path):
        subgraph = -1
        line = LINE.strip("\n").split("\t")
        if line[0] == "S":
            idx = line[1]
            subgraph = vertex_to_subgraph[idx]
        elif line[0] == "L":
            idx1, idx2 = line[1], line[3]
            if vertex_to_subgraph[idx1] == vertex_to_subgraph[idx2]:
                subgraph = vertex_to_subgraph[idx1]
        elif line[0] == "P":
            idx = line[2].split(",")[0][:-1]
            subgraph = vertex_to_subgraph[idx]
            # TODO remove this for
            for v in line[2].split(",")[1:]:
                idx = v[:-1]
                assert subgraph == vertex_to_subgraph[idx]
        if subgraph == -1:
            continue
        ofiles[subgraph].write(LINE)
    for ofile in ofiles:
        ofile.close()


if __name__ == "__main__":
    main()
