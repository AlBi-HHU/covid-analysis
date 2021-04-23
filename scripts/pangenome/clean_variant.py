import vcfpy
import json
import math

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict, OrderedDict
import statistics
import sys


sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import read_node2len, strand_bias


def main(
    in_vcf,
    supportFile,
    node2len_path,
    min_cov,
    rvt,
    out_vcf,
):

    support = json.load(open(supportFile, "r"))
    nodeSupport = support["nodes"]

    edgeSupport = defaultdict(lambda: defaultdict(int))
    for key, value in support["edges"].items():
        edgeSupport[frozenset(key.split("_"))] = value

    node2len = read_node2len(node2len_path)

    reader = vcfpy.Reader.from_path(in_vcf)

    pos2var = defaultdict(list)
    for record in reader:
        vcoverage = record.INFO["VCOVT"]
        pos2var[(record.POS, tuple(record.REF), tuple(record.ALT))].append(
            (vcoverage, record, record.INFO["VCOVT"])
        )

    # Add info line
    header = create_header(reader.header)

    # Add filter line
    add_header_filter(header)

    variantsToWrite = []

    for key, values in pos2var.items():
        variant = values[0][1]

        if len(values) != 1:
            values.sort(key=lambda x: x[0], reverse=True)
            vcoverage = values[0][0]
            variant = values[0][1]

        vsup_f, vsup_r = compute_support(
            variant.INFO["VARPATH"], node2len, nodeSupport, edgeSupport, strict=True
        )
        rsup_f, rsup_r = compute_support(
            variant.INFO["REFPATH"], node2len, nodeSupport, edgeSupport, strict=True
        )
        vsup = vsup_f + vsup_r
        rsup = rsup_f + rsup_r

        if (min(rsup_f, rsup_r) + min(vsup_f, vsup_r)) == 0:
            variant.INFO["CORHETRATIO"] = -1.0
        else:
            cor_vsup = min(vsup_f, vsup_r) * 2
            cor_rsup = min(rsup_f, rsup_r) * 2
            variant.INFO["CORHETRATIO"] = cor_vsup / (cor_vsup + cor_rsup)

        if vsup_f != float("nan") and vsup_r != float("nan"):
            vcoverage = vsup_f + vsup_r
        if rsup_f != float("nan") and rsup_r != float("nan"):
            rcoverage = rsup_f + rsup_r

        filters = list()
        if vcoverage < min_cov:
            filters.append("VCoverage")
        elif not math.isnan(vsup):
            if strand_bias(variant, component="V", sb_max_threshold=float(snakemake.config["pangenomeMaxSB"])):
                filters.append("VStrandBias")

        if rcoverage < min_cov:
            filters.append("RCoverage")
        elif not math.isnan(rsup):
            if strand_bias(variant, component="R", sb_max_threshold=float(snakemake.config["pangenomeMaxSB"])):
                filters.append("RStrandBias")

        if (vsup / (vsup + rsup)) < rvt:
            filters.append("LowRVT")

        if len(filters) == 0:
            filters.append("PASS")

        variant.FILTER = filters

        variant = rebind_info(variant)

        variantsToWrite.append(variant)

    # drop additional fields
    header = vcfpy.header_without_lines(
        header,
        [
            ("INFO", "RCOVT"),
            ("INFO", "VCOVT"),
            ("INFO", "RCOVF"),
            ("INFO", "RCOVR"),
            ("INFO", "VCOVF"),
            ("INFO", "VCOVR"),
            ("INFO", "OLD_VARIANT"),
            ("INFO", "OLD_CLUMPED"),
            ("INFO", "VARPATH"),
            ("INFO", "REFPATH"),
        ],
    )

    writer = vcfpy.Writer.from_path(out_vcf, header)
    for variant in variantsToWrite:
        writer.write_record(variant)


def compute_support(nodes, node2len, node_support, edge_support, strict=False):
    nodes = nodes.split("_")

    if len(nodes) <= 2:
        forward = (
            edge_support[frozenset(nodes)]["forward"]
            if "forward" in edge_support[frozenset(nodes)]
            else 0
        )
        reverse = (
            edge_support[frozenset(nodes)]["reverse"]
            if "reverse" in edge_support[frozenset(nodes)]
            else 0
        )
        return (forward, reverse)

    all_supports_f = 0
    all_supports_r = 0
    path_len = 0

    key = "strict" if strict else "lenient"

    for node in nodes[1:-1]:
        # Take maximum of coverage across all positions
        all_supports_f += statistics.mean(
            node_support[node]["forward"][pos][key]
            for pos in range(len(node_support[node]["forward"]))
        )
        all_supports_r += statistics.mean(
            node_support[node]["reverse"][pos][key]
            for pos in range(len(node_support[node]["reverse"]))
        )
        path_len += node2len[node]

    return (all_supports_f / path_len, all_supports_r / path_len)


def rebind_info(record):
    old_info = record.INFO
    new_info = OrderedDict()

    new_info["BUBBLEID"] = str(old_info["BUBBLEID"])
    new_info["RPATH"] = old_info["REFPATH"]
    new_info["VPATH"] = old_info["VARPATH"]

    new_info["RCOV"] = [old_info["RCOVF"], old_info["RCOVR"]]
    new_info["VCOV"] = [old_info["VCOVF"], old_info["VCOVR"]]

    new_info["CORHETRATIO"] = old_info["CORHETRATIO"]

    record.INFO = new_info

    return record


def create_header(header):

    header.add_info_line(
        {
            "ID": "RPATH",
            "Type": "String",
            "Number": "1",
            "Description": "Ids of nodes that make up the reference path (separate by _)",
        }
    )

    header.add_info_line(
        {
            "ID": "VCOV",
            "Type": "Float",
            "Number": "2",
            "Description": "Coverage for the Variant Allele (Tuple Forward,Backward)",
        }
    )

    header.add_info_line(
        {
            "ID": "RCOV",
            "Type": "Float",
            "Number": "2",
            "Description": "Coverage for the Reference Allele (Tuple Forward,Backward)",
        }
    )

    header.add_info_line(
        {
            "ID": "VPATH",
            "Type": "String",
            "Number": "1",
            "Description": "Ids of nodes that make up the variant path (separate by _)",
        }
    )

    header.add_info_line(
        {
            "ID": "CORHETRATIO",
            "Type": "Float",
            "Number": "1",
            "Description": "Corrected heterozygote ratio variant support divide by variant plus reference support",
        }
    )

    return header


def add_header_filter(header):
    header.add_filter_line(
        {
            "ID": "VCoverage",
            "Description": "Variant coverage is lower than minimal coverage",
        }
    )

    header.add_filter_line(
        {
            "ID": "VStrandBias",
            "Description": "We notice a strand bias in coverage of this variant",
        }
    )

    header.add_filter_line(
        {
            "ID": "RCoverage",
            "Description": "Reference coverage is lower than minimal coverage",
        }
    )

    header.add_filter_line(
        {
            "ID": "RStrandBias",
            "Description": "We notice a strand bias in coverage of the reference emissions",
        }
    )

    header.add_filter_line(
        {
            "ID": "LowRVT",
            "Description": "RVT is below threshold",
        }
    )


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.input["support"],
        snakemake.input["node2len"],
        snakemake.config["pangenomeVarMinCov"],
        snakemake.config["pangenomeRVTTSupport"],
        snakemake.output[0],
    )
else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
