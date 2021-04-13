import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import get_node2seq
import re
import json

from collections import defaultdict


def main(alignment, pangenome, output):

    # get the sequence associate to node
    node2seq = get_node2seq(pangenome)

    node2base_cov = {}

    for n_id, seq in node2seq.items():
        node2base_cov[n_id] = {"forward": [], "reverse": []}
        for pos in range(len(seq)):
            node2base_cov[n_id]["forward"].append({"strict": 0, "lenient": 0})
            node2base_cov[n_id]["reverse"].append({"strict": 0, "lenient": 0})

    edge2cov = defaultdict(lambda: defaultdict(int))

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
                data[5][1:].split(">")
                if orientation == "forward"
                else data[5][1:].split("<")[::-1]
            )

            # 7 -> Path Start, 8 -> Path End, 6 -> Path Length
            path_start = (
                int(data[7])
                if (orientation == "forward")
                else int(data[6]) - int(data[8])
            )

            cigar = data[-1].split(":")[2]
            cigar = re.findall("(\d+[MNDI])", cigar)

            logFile.write(",".join(cigar) + "\n")

            # If the last instruction is an I, discard it, it has no information
            if cigar[-1][-1] == "I":
                cigar = cigar[:-1]

            current_node = path.pop(0)
            current_node_len = len(node2seq[current_node])
            position_on_node = (
                path_start
                if orientation == "forward"
                else ((current_node_len - path_start) - 1)
            )
            logFile.write(
                "Starting on node {} at position {}".format(
                    current_node, position_on_node
                )
                + "\n"
            )

            # Memoize previous CIGAR instructions for stricter values
            lastInstructions = ["None", "None", "None"]
            lastPositions = [-1, -1, -1]
            lastNodes = [-1, -1, -1]

            for (*length, instruction) in cigar:
                length = int("".join(length))

                for posInInstruction in range(length):
                    if orientation == "forward":
                        if position_on_node >= current_node_len:
                            if len(path) == 0:
                                logFile.write("Warning, path is exhausted!" + "\n")
                                logFile.write("ReadID: {}".format(data[0]) + "\n")
                                exit(59)
                            current_node = path.pop(0)
                            current_node_len = len(node2seq[current_node])
                            position_on_node = 0
                    else:
                        if position_on_node < 0:
                            if len(path) == 0:
                                logFile.write("Warning, path is exhausted!" + "\n")
                                logFile.write("ReadID: {}".format(data[0]) + "\n")
                                exit(59)
                            current_node = path.pop(0)
                            current_node_len = len(node2seq[current_node])
                            position_on_node = current_node_len - 1
                            logFile.write("{},{}".format(length, instruction) + "\n")
                            logFile.write(",".join(path) + "\n")

                    # print(current_node, len(node2base_cov[current_node][orientation]), position_on_node, posInInstruction, instruction, length, cigar, line)

                    if instruction == "M":
                        node2base_cov[current_node][orientation][position_on_node][
                            "lenient"
                        ] += 1

                    lastInstructions.pop(-1)
                    lastPositions.pop(-1)
                    lastNodes.pop(-1)

                    lastInstructions.insert(0, instruction)
                    lastPositions.insert(0, position_on_node)
                    lastNodes.insert(0, current_node)

                    # if we have three matches, the middle match is REALLY supported
                    if (
                        lastInstructions[0]
                        == lastInstructions[1]
                        == lastInstructions[2]
                        == "M"
                    ):
                        node2base_cov[lastNodes[1]][orientation][lastPositions[1]][
                            "strict"
                        ] += 1

                    # test if we are at end of a node
                    if (orientation == "forward" and position_on_node == 0) or (
                        orientation == "reverse"
                        and position_on_node == current_node_len - 1
                    ):
                        if lastInstructions[0] == lastInstructions[1] == "M":
                            normalize_key = (
                                (lastNodes[0], lastNodes[1])
                                if lastNodes[0] < lastNodes[1]
                                else (lastNodes[1], lastNodes[0])
                            )
                            edge2cov["_".join(normalize_key)][orientation] += 1

                    if instruction != "I":
                        if orientation == "forward":
                            position_on_node += 1
                        else:
                            position_on_node -= 1

        json.dump({"nodes": node2base_cov, "edges": edge2cov}, outFile)


if "snakemake" in locals():
    main(
        snakemake.input["alignment"], snakemake.input["pangenome"], snakemake.output[0]
    )

else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3])
