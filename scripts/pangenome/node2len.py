import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import get_node2seq


def main(gfa, output):
    node2seq = get_node2seq(gfa)

    with open(output, "w") as out:
        print(f"node,length", file=out)

        for (node, seq) in node2seq.items():
            print(f"{node},{len(seq)}", file=out)


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.output["data"],
    )
else:
    main(sys.argv[1], sys.argv[2])
