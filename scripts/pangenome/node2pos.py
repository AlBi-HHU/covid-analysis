import re
import csv

from collections import defaultdict

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure


from shared import get_node2seq


def main(graph_path, reference_mapping, node2pos_path):

    ref_path = list(get_ref_path(reference_mapping))

    # Get sequence of each node
    node2seq = get_node2seq(graph_path)

    # Iterate over path to get position by cumulative length
    node2ref_pos = dict()
    node2ori = dict()
    cumu_len = 0
    for node in ref_path:
        node2ori[node[1:]] = node[0]
        node2ref_pos[node[1:]] = cumu_len
        cumu_len += len(node2seq[node[1:]])

    # Write result
    with open(node2pos_path, "w") as csv_fh:
        print("node,pos,ori", file=csv_fh)
        for node, pos in sorted(node2ref_pos.items(), key=lambda x: x[1]):
            print("{},{},{}".format(node, pos, node2ori[node]), file=csv_fh)


def get_ref_path(path):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            nodes = re.findall(r"([<>][^<>]+)", row[5])

            for n in nodes:
                yield n



if "snakemake" in locals():
    main(
        snakemake.input["graph"],
        snakemake.input["mappings"],
        snakemake.output["node_pos_on_ref"],
    )
else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3])
