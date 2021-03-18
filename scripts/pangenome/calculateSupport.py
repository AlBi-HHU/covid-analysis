import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import get_node2seq
import re
import json


def main(alignment, pangenome, output):


    # get the sequence associate to node
    node2seq = get_node2seq(pangenome)

    node2base_cov = {n_id: {"forward": [0] * len(seq), "reverse": [0] * len(seq)} for n_id, seq in node2seq.items()}

    with open(alignment, "r") as alignmentFile, open(
        snakemake.log[0], "w"
    ) as logFile, open(output, "w") as outFile:
        # For each line in the alignment file
        for line in alignmentFile:
            # We split based on whitespace
            data = line.split()

            # Split Read, the first Symbol in path is only used as an indicator for the orientation
            orientation = "forward" if data[5][0] == ">" else "reverse"
            path = list(
                data[5][1:].split(">")[::-1]
                if orientation == "forward"
                else data[5][1:].split("<")
            )

            path_length = sum(len(node2seq[n_id]) for n_id in path)
            path_start = (
                int(data[7]) if (orientation == 'forward') else path_length - int(data[8])
            )

            cigar = data[-1].split(":")[2]
            cigar = re.findall("(\d+[MNDI])", cigar)

            # If the last instruction is an I, discard it, it has no information
            if (cigar[-1][-1] == "I"):
                cigar = cigar[:-1]

            current_node = path.pop(0)
            current_node_len = len(node2seq[current_node])
            position_on_node = path_start

            for (*length, instruction) in cigar:

                length = int("".join(length))

                for _ in range(length):
                    if position_on_node >= current_node_len:
                        if len(path) == 0:
                            print("Warning, path is exhausted!")
                            exit(59)
                        current_node = path.pop(0)
                        current_node_len = len(node2seq[current_node])
                        position_on_node = 0

                    if instruction == "M":
                        node2base_cov[current_node][orientation][position_on_node] += 1

                    if instruction != "I":
                        position_on_node += 1

        json.dump(node2base_cov, outFile)


if "snakemake" in locals():
    main(
        snakemake.input["alignment"], snakemake.input["pangenome"], snakemake.output[0]
    )

else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3])
