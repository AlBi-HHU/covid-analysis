import vcfpy
import json
import math

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict, OrderedDict
import statistics
import sys

from scipy.stats import binom

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import read_node2len


def main(
    in_vcf,
    supportFile,
    node2len_path,
    min_cov,
    rvt,
    th_sbiais,
    out_vcf,
):

    nodeSupport = json.load(open(supportFile, "r"))
    node2len = read_node2len(node2len_path)

    reader = vcfpy.Reader.from_path(in_vcf)

    header = clean_header_info(reader.header)

    pos2var = defaultdict(list)
    for record in reader:
        coverage = record.INFO["VCOV"]
        pos2var[(record.POS, tuple(record.REF), tuple(record.ALT))].append(
            (coverage, record, record.INFO["VCOV"])
        )

    # Add info line
    add_header_info(header)

    # Add filter line
    add_header_filter(header)

    writer = vcfpy.Writer.from_path(out_vcf, reader.header)

    for key, values in pos2var.items():
        variant = values[0][1]

        if len(values) != 1:
            values.sort(key=lambda x: x[0], reverse=True)
            coverage = values[0][0]
            variant = values[0][1]

        vsup_f, vsup_r = compute_support(variant.INFO["VARPATH"], node2len, nodeSupport,strict=True)
        rsup_f, rsup_r = compute_support(variant.INFO["REFPATH"], node2len, nodeSupport,strict=True)
        vsup = vsup_f + vsup_r
        rsup = rsup_f + rsup_r

        variant.INFO["VSUP"] = vsup
        variant.INFO["RSUP"] = rsup
        variant.INFO["VSUPF"] = vsup_f
        variant.INFO["VSUPR"] = vsup_r
        variant.INFO["RSUPF"] = rsup_f
        variant.INFO["RSUPR"] = rsup_r


        if (min(rsup_f, rsup_r) + min(vsup_f, vsup_r)) == 0:
            variant.INFO["CORHETRATIO"] = -1.0
        else:
            cor_vsup = min(vsup_f, vsup_r) * 2
            cor_rsup = min(rsup_f, rsup_r) * 2
            variant.INFO["CORHETRATIO"] = cor_vsup / (cor_vsup + cor_rsup)

        if vsup_f != float("nan") and rsup_f != float("nan"):
            coverage = vsup_f + vsup_r

        filters = list()
        if coverage < min_cov:
            filters.append("Coverage")
        elif not math.isnan(vsup):
            if strand_bias(variant, th_sbiais, "SUP"):
                filters.append("StrandBias")

            if (vsup + rsup) == 0 or (vsup / (vsup + rsup)) < rvt:
                filters.append("NoRealSupport")

        if len(filters) == 0:
            filters.append("PASS")

        variant.FILTER = filters

        variant = rebind_info(variant)

        writer.write_record(variant)


def clean_header_info(header):
    del_position = reversed([p for (p, line) in enumerate(header.lines) if isinstance(line, vcfpy.InfoHeaderLine)])

    for pos in del_position:
        del header.lines[pos]

    return header


def compute_support(nodes, node2len, node_support,strict=False):
    nodes = nodes.split("_")

    if len(nodes) <= 2:
        return (float("nan"), float("nan"))

    all_supports_f = 0
    all_supports_r = 0
    path_len = 0

    key = 'strict' if strict else 'lenient'

    for node in nodes[1:-1]:
        #Take maximum of coverage across all positions
        all_supports_f += statistics.mean(node_support[node]["forward"][pos][key] for pos in range(len(node_support[node]["forward"])))
        all_supports_r += statistics.mean(node_support[node]["reverse"][pos][key] for pos in range(len(node_support[node]["reverse"])))
        path_len += node2len[node]

    return (all_supports_f / path_len, all_supports_r / path_len)


def strand_bias(record, th_sbias, text="COV"):
    vcov = float(record.INFO[f"V{text}"])
    vcov_forward = float(record.INFO[f"V{text}F"])
    vcov_reverse = float(record.INFO[f"V{text}R"])

    ratio = min(vcov_forward, vcov_reverse) / vcov
    if ratio > th_sbias:
        return False
    else:
        return True


def rebind_info(record):
    old_info = record.INFO
    new_info = OrderedDict()

    new_info["BUBBLEID"] = str(old_info["BUBBLEID"])
    new_info["RPATH"] = old_info["REFPATH"]
    new_info["VPATH"] = old_info["VARPATH"]

    new_info["RCOV"] = [old_info["RCOVF"], old_info["RCOVR"]]
    new_info["VCOV"] = [old_info["VCOVF"], old_info["VCOVR"]]

    new_info["RSUP"] = [old_info["RSUPF"], old_info["RSUPR"]]
    new_info["VSUP"] = [old_info["VSUPF"], old_info["VSUPR"]]

    new_info["CORHETRATIO"] = old_info["CORHETRATIO"]

    record.INFO = new_info

    return record


def add_header_info(header):
    header.add_info_line(
        {
            "ID": "BUBBLEID",
            "Type": "String",
            "Number": "1",
            "Description": "Id of bubble in pangenome",
        }
    )

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
            "ID": "VPATH",
            "Type": "String",
            "Number": "1",
            "Description": "Ids of nodes that make up the variant path (separate by _)",
        }
    )

    header.add_info_line(
        {
            "ID": "RCOV",
            "Type": "Float",
            "Number": "2",
            "Description": "Reference coverage on each strand",
        }
    )

    header.add_info_line(
        {
            "ID": "VCOV",
            "Type": "Float",
            "Number": "2",
            "Description": "Variant support on each strand",
        }
    )

    header.add_info_line(
        {
            "ID": "RSUP",
            "Type": "Float",
            "Number": "2",
            "Description": "Reference support on each strand",
        }
    )

    header.add_info_line(
        {
            "ID": "VSUP",
            "Type": "Float",
            "Number": "2",
            "Description": "Variant support on each strand",
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


def add_header_filter(header):
    header.add_filter_line(
        {
            "ID": "Coverage",
            "Description": "Total coverage is lower than minimal coverage",
        }
    )

    header.add_filter_line(
        {
            "ID": "StrandBias",
            "Description": "We notice a strand bias in coverage of this variant",
        }
    )

    header.add_filter_line(
        {
            "ID": "NoRealSupport",
            "Description": "This variant isn't realy support",
        }
    )


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.input["support"],
        snakemake.input["node2len"],
        snakemake.config["pangenomeVarMinCov"],
        snakemake.config["pangenomeRVTTSupport"],
        float(snakemake.config["pangenomeStrandBiais"]),
        snakemake.output[0],
    )
else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
