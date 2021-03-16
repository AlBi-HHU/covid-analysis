import vcfpy
import json

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict

import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import read_node2len


def main(in_vcf, supportFile, node2len_path, min_cov, rvt, out_vcf):

    nodeSupport = json.load(open(supportFile, "r"))
    node2len = read_node2len(node2len_path)

    reader = vcfpy.Reader.from_path(in_vcf)

    header = reader.header

    pos2var = defaultdict(list)
    for record in reader:
        coverage = record.INFO["VCOV"] + record.INFO["RCOV"]
        pos2var[(record.POS, tuple(record.REF), tuple(record.ALT))].append(
            (coverage, record, record.INFO["VCOV"])
        )

    header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "LRS"),
                (
                    "Description",
                    "Real support below {} percent".format(
                        snakemake.config["pagenomeCutoffRealSupport"]
                    ),
                ),
            ]
        )
    )
    header.add_info_line(
        {
            "ID": "MULTIPLE",
            "Type": "Flag",
            "Number": "1",
            "Description": "Pangenome found multiple variant at this position",
        }
    )
    header.add_info_line(
        {
            "ID": "VSUP",
            "Type": "String",
            "Number": "1",
            "Description": "Variant support",
        }
    )
    header.add_info_line(
        {
            "ID": "RSUP",
            "Type": "String",
            "Number": "1",
            "Description": "Reference support",
        }
    )
    writer = vcfpy.Writer.from_path(out_vcf, reader.header)

    for key, values in pos2var.items():

        variant = values[0][1]

        if len(values) != 1:
            values.sort(key=lambda x: x[0], reverse=True)
            coverage = values[0][0]
            variant = values[0][1]
            variant.INFO["MULTIPLE"] = True

        vsup = compute_support(variant.INFO["VARPATH"])
        rsup = compute_support(variant.INFO["REFPATH"])

        coverage = vsup + rsup

        if coverage > min_cov and (vsup / (vsup + rsup)) > rvt:
            writer.write_record(variant)


def compute_support(nodes, node2len, node_support):
    nodes = nodes.spilt(",")
    all_supports = 0
    path_len = 0

    for node in nodes:
        all_supports += node_support[node][1]
        path_len += node2len[node]

    return all_supports / path_len


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.input["support"],
        snakemake.input["node2len"],
        snakemake.config["pangenomeVarMinCov"],
        snakemake.config["pangenomeRVTThreshold"],
        snakemake.output[0],
    )
else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
