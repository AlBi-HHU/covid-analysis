import vcfpy
import json
import math

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict

import sys

from scipy.stats import binom

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import read_node2len


def main(in_vcf, supportFile, node2len_path, min_cov, rvt, th_sbiais, th_sb_cov, th_sb_pval, out_vcf):

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

    # Add info line
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
            "ID": "VSUPF",
            "Type": "String",
            "Number": "1",
            "Description": "Variant support forward",
        }
    )
    header.add_info_line(
        {
            "ID": "VSUPR",
            "Type": "String",
            "Number": "1",
            "Description": "Variant support reverse",
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
    header.add_info_line(
        {
            "ID": "RSUPF",
            "Type": "String",
            "Number": "1",
            "Description": "Reference support forward",
        }
    )
    header.add_info_line(
        {
            "ID": "RSUPR",
            "Type": "String",
            "Number": "1",
            "Description": "Reference support reverse",
        }
    )

    # Add filter line
    header.add_filter_line(
        {
            "ID": "Coverage",
            "Description": "Total coverage is lower than minimal coverage",
        }
    )

    header.add_filter_line(
        {
            "ID": "StrandBiais",
            "Description": "We notice a strand biais in coverage of this variant",
        }
    )

    header.add_filter_line(
        {
            "ID": "NoRealSupport",
            "Description": "This variant isn't realy support",
        }
    )

    header.add_filter_line(
        {
            "ID": "StrandBiaisRealSupport",
            "Description": "We notice a strand biais in coverage of this variant",
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

        vsup_f, vsup_r = compute_support(variant.INFO["VARPATH"], node2len, nodeSupport)
        rsup_f, rsup_r = compute_support(variant.INFO["REFPATH"], node2len, nodeSupport)
        vsup = vsup_f + vsup_r
        rsup = rsup_f + rsup_r

        variant.INFO["VSUP"] = vsup
        variant.INFO["RSUP"] = rsup
        variant.INFO["VSUPF"] = vsup_f
        variant.INFO["VSUPR"] = vsup_r
        variant.INFO["RSUPF"] = rsup_f
        variant.INFO["RSUPR"] = rsup_r

        if vsup_f != float("nan") and rsup_f != float("nan"):
            coverage = vsup_f + vsup_r + rsup_f + rsup_r

        filters = list()
        if coverage < min_cov:
            filters.append("Coverage")

        if strand_biais(variant, th_sbiais, th_sb_cov, th_sb_pval, "COV"):
            filters.append("StrandBiais")

        if not math.isnan(vsup):
            if strand_biais(variant, th_sbiais, th_sb_cov, th_sb_pval, "SUP"):
                filters.append("StrandBiaisRealSupport")

            if (vsup / (vsup + rsup)) < rvt:
                filters.append("NoRealSupport")

        if len(filters) == 0:
            filters.append("PASS")

        variant.FILTER = filters

        writer.write_record(variant)


def compute_support(nodes, node2len, node_support):
    nodes = nodes.split("_")

    if len(nodes) <= 2:
        return (float("nan"), float("nan"))

    all_supports_f = 0
    all_supports_r = 0
    path_len = 0

    for node in nodes[1:-1]:
        all_supports_f += max(node_support[node]["forward"])
        all_supports_r += max(node_support[node]["reverse"])
        path_len += node2len[node]

    return (all_supports_f / path_len, all_supports_r / path_len)


def strand_biais(record, th_sbiais, th_sb_cov, th_sb_pval, text="COV"):
    rcov = float(record.INFO[f"R{text}"])
    vcov = float(record.INFO[f"V{text}"])
    vcov_forward = float(record.INFO[f"V{text}F"])
    vcov_reverse = float(record.INFO[f"V{text}R"])
    tcov = rcov + vcov

    if vcov < th_sb_cov:
        pval = binom.pmf(vcov_forward, vcov, 0.5)
        if pval > th_sb_pval:
            return False
        else:
            return True
    else:
        ratio = min(vcov_forward, vcov_reverse) / tcov
        if ratio > th_sbiais:
            return False
        else:
            return True


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.input["support"],
        snakemake.input["node2len"],
        snakemake.config["pangenomeVarMinCov"],
        snakemake.config["pangenomeRVTTSupport"],
        float(snakemake.config["pangenomeStrandBiais"]),
        int(snakemake.config['pangenomePValCoverage']),
        float(snakemake.config['pangenomePValCutoff']),
        snakemake.output[0],
    )
else:
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
