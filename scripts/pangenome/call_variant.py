import csv
import json
import networkx
import sys
from collections import defaultdict, Counter

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import *  # TODO: import only required modules


def main(
    pangenome_path,
    bubble_path,
    reads_mapping,
    node2pos_path,
    rvt_threshold,
    min_cov_factor,
    min_ends_cov,
    output_path,
):

    # Get ref information
    node2pos = dict()
    node2ori = dict()
    ref_path = list()

    # Read node position on reference
    with open(node2pos_path) as fh:
        reader = csv.reader(fh, delimiter=",")
        next(reader)
        for row in reader:
            node2pos[row[0]] = int(row[1])
            node2ori[row[0]] = row[2]
            ref_path.append(row[0])

    ref_nodes = set(ref_path)
    ref_edge = {(ref_path[i], ref_path[i + 1]) for i in range(len(ref_path) - 1)}

    # get the sequence associate to node
    node2seq = get_node2seq(pangenome_path)
    graph = gfa2networkx(pangenome_path)

    # Get node coverage
    node2base = defaultdict(lambda: Counter())
    edge2cov = Counter()
    paths = defaultdict(list)

    parse_gaf(reads_mapping, paths, node2base, edge2cov, node2seq)

    node2cov = defaultdict(lambda: Counter())
    for (node, base) in node2base.items():
        node2cov[node][True] = base[True] / len(node2seq[node])
        node2cov[node][False] = base[False] / len(node2seq[node])
        node2cov[node]["all"] = (base[True] + base[False]) / len(node2seq[node])

    # Read bubble
    simple_bubble = set()
    bubbles = dict()
    chains = json.load(open(bubble_path))
    for (chain_id, chain) in chains.items():
        if "parent_sb" in chain:
            super_bubble = chain["parent_sb"]
        else:
            super_bubble = None

        for bubble in chain["bubbles"]:
            id = bubble["id"]
            del bubble["id"]

            if super_bubble is not None:
                bubble["parent"] = super_bubble

            if bubble["type"] == "simple" or bubble["type"] == "insertion":
                simple_bubble.add(id)
            bubbles[id] = bubble

    variants = list()
    for b_id in bubbles:

        # Get bubble with ends in ref
        skip_bubble = False
        bubble = bubbles[b_id]
        while not all((end in ref_nodes for end in bubble["ends"])):
            if "parent" in bubble:
                bubble = bubbles[bubble["parent"]]
            else:
                skip_bubble = True
                break

        if skip_bubble:
            continue

        # Compute coverage around bubble
        ends_cov = min((node2cov[node]["all"] for node in bubble["ends"]))

        # If bubble have ref node with low coverage it's a variant bubble
        if ends_cov < min_ends_cov:
            continue

        subgraph = graph.subgraph(bubble["ends"] + bubble["inside"]).copy()

        # Annotate node
        not_cov_nodes = set()
        not_cov_edges = set()
        for (node, data) in subgraph.nodes.items():
            if node2cov[node]["all"] < ends_cov * min_cov_factor:
                not_cov_nodes.add(node)

        # Annotate edge
        for (edge, data) in subgraph.edges.items():
            if edge2cov[frozenset(edge)] < ends_cov * min_cov_factor:
                not_cov_edges.add(edge)

        # Found reference path in bubble
        ref_path = get_paths(
            subgraph.subgraph((n for n in subgraph.nodes() if n in ref_nodes)),
            bubble["ends"],
        )
        if len(ref_path) == 0:
            continue
        elif len(ref_path) == 1:
            ref_path = ref_path[0]
        else:
            ref_path = sorted(ref_path, key=lambda x: len(x), reverse=True)[0]

        ref_seq = "".join(node2seq[node] for node in ref_path)
        (ref_cov, ref_covf, ref_covr) = path_coverage(
            ref_path, edge2cov, node2cov, node2seq
        )

        # Clean not covered node and edge
        subgraph.remove_nodes_from(not_cov_nodes)
        subgraph.remove_edges_from(not_cov_edges)

        # Found variant paths
        var_paths = list()
        for path in get_paths(subgraph, bubble["ends"]):
            prev = None
            for node in path:
                if node not in ref_nodes:
                    var_paths.append(path)
                    break

                if prev is not None and (prev, node) not in ref_edge:
                    var_paths.append(path)
                    break

                prev = node

        for var_path in var_paths:
            (var_cov, var_covf, var_covr) = path_coverage(
                var_path, edge2cov, node2cov, node2seq
            )

            var_seq = "".join(node2seq[node] for node in var_path)

            variants.append(
                (
                    ref_seq,
                    var_seq,
                    node2pos[ref_path[0]],
                    ref_cov,
                    ref_covf,
                    ref_covr,
                    var_cov,
                    var_covf,
                    var_covr,
                    b_id,
                    ref_path[1:-1],
                    var_path[1:-1],
                )
            )

    # set data required by vcf format
    ref_name = "MN908947.3"
    ref_length = 29903

    # Write variant
    with open(output_path, "w") as fh:
        vcf_header(fh, ref_name, ref_length)

        for (
            r_seq,
            v_seq,
            pos,
            ref_cov,
            ref_covf,
            ref_covr,
            var_cov,
            var_covf,
            var_covr,
            bubble_id,
            ref_path,
            var_path,
        ) in variants:
            rvt = var_cov / (ref_cov + var_cov)
            if rvt >= rvt_threshold:
                print(
                    "{}\t{}\t.\t{}\t{}\t.\t.\tRCOV={:.4f};RCOVF={:.4f};RCOVR={:.4f};VCOV={:.4f};VCOVF={:.4f};VCOVR={:.4f};BUBBLEID={};REALPATH={};VARPATH={}".format(
                        ref_name,
                        pos + 1,
                        r_seq,
                        v_seq,
                        ref_cov,
                        ref_covf,
                        ref_covr,
                        var_cov,
                        var_covf,
                        var_covr,
                        bubble_id,
                        ",".join(ref_path),
                        ",".join(var_path),
                    ),
                    file=fh,
                )


def path_coverage(path, edge2cov, node2cov, node2seq):
    if len(path) < 2:
        return (0, 0, 0)
    elif len(path) == 2:
        if frozenset(path) in edge2cov:
            return (
                edge2cov[frozenset(path)],
                edge2cov[frozenset(path)],
                edge2cov[frozenset(path)],
            )
        else:
            return (0, 0, 0)

    path_len = sum(len(node2seq[x]) for x in path[1:-1])

    # Don't take ends whne compute var_cov
    return (
        sum(node2cov[node]["all"] for node in path[1:-1]) / path_len,
        sum(node2cov[node][True] for node in path[1:-1]) / path_len,
        sum(node2cov[node][False] for node in path[1:-1]) / path_len,
    )


def get_paths(graph, ends):
    paths = list(networkx.all_simple_paths(graph, ends[0], ends[1]))
    if len(paths) == 0:
        paths = list(networkx.all_simple_paths(graph, ends[1], ends[0]))

    return paths


def vcf_header(fh, ref_name, length):
    print("##fileformat=VCFv4.2", file=fh)
    print(
        '##INFO=<ID=RCOV,Number=1,Type=Float,Description="Coverage of reference path">',
        file=fh,
    )
    print(
        '##INFO=<ID=RCOVF,Number=1,Type=Float,Description="Coverage of reference path forward strand">',
        file=fh,
    )
    print(
        '##INFO=<ID=RCOVR,Number=1,Type=Float,Description="Coverage of reference path reverse strand ">',
        file=fh,
    )
    print(
        '##INFO=<ID=VCOV,Number=1,Type=Float,Description="Coverage of variant path">',
        file=fh,
    )
    print(
        '##INFO=<ID=VCOVF,Number=1,Type=Float,Description="Coverage of variant path forward strand">',
        file=fh,
    )
    print(
        '##INFO=<ID=VCOVR,Number=1,Type=Float,Description="Coverage of variant path reverse strand">',
        file=fh,
    )
    print(
        '##INFO=<ID=BUBBLEID,Number=1,Type=Integer,Description="Id of bubble in pangenome">',
        file=fh,
    )
    print(
        '##INFO=<ID=REFPATH,Number=1,Type=String,Description="Ids of nodes that make up the reference path">',
        file=fh,
    )
    print(
        '##INFO=<ID=VARPATH,Number=1,Type=String,Description="Ids of nodes that make up the variant path">',
        file=fh,
    )
    print("##contig=<ID={},length={}>".format(ref_name, length), file=fh)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=fh)


def sequence_from_node(oris, nodes, node2seq, previous_ori):
    seq = ""
    for ori, node in zip(oris, nodes):
        if ori == ">":
            seq += node2seq[node]
        else:
            seq += rev_comp(node2seq[node])

    return seq


def pos_of_diff(path, reference):
    begin = path[0]

    try:
        index_ref = reference.index(begin)
    except ValueError:
        index_ref = None

    ref_set = set(reference)
    for i_not_ref in (
        index for (index, node) in enumerate(path) if node not in set(reference)
    ):
        begin_break = None
        # found reference node before variant
        for i in range(i_not_ref, 0, -1):
            if path[i] in ref_set:
                begin_break = i
                break

        # found reference node after variant
        end_break = None
        for i in range(i_not_ref, len(path)):
            if path[i] in ref_set:
                end_break = i
                break

        if index_ref is None or begin_break is None or end_break is None:
            yield None
        else:
            yield (begin_break, end_break + 1)


def gfa2networkx(gfa_path):
    graph = networkx.DiGraph()

    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("L"):
                line = line.split("\t")
                graph.add_edge(line[1], line[3])

    return graph


if "snakemake" in locals():
    rvt = float(snakemake.config["pangenomeRVTTPath"])
    min_cov_factor = float(snakemake.config["pangenomeMinCovFactor"])
    min_ends_cov = int(snakemake.config["pangenomeMinEndsCovBubble"])
    main(
        snakemake.input["pangenome"],
        snakemake.input["bubble"],
        snakemake.input["reads"],
        snakemake.input["node2pos"],
        rvt,
        min_cov_factor,
        min_ends_cov,
        snakemake.output["variant"],
    )
else:
    import sys

    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        float(sys.argv[4]),
        float(sys.argv[5]),
        sys.argv[6],
    )
