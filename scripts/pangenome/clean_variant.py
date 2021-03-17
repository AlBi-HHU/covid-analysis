import vcfpy
import json

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
            "ID": "RSUP",
            "Type": "String",
            "Number": "1",
            "Description": "Reference support",
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

    writer = vcfpy.Writer.from_path(out_vcf, reader.header)

    for key, values in pos2var.items():
        variant = values[0][1]

        if len(values) != 1:
            values.sort(key=lambda x: x[0], reverse=True)
            coverage = values[0][0]
            variant = values[0][1]
            variant.INFO["MULTIPLE"] = True

        vsup = compute_support(variant.INFO["VARPATH"], node2len, nodeSupport)
        rsup = compute_support(variant.INFO["REFPATH"], node2len, nodeSupport)

        variant.INFO["VSUP"] = vsup
        variant.INFO["RSUP"] = rsup

        if vsup != float("nan") and rsup != float("nan"):
            coverage = vsup + rsup

        filters = list()
        if coverage < min_cov:
            filters.append("Coverage")

        if strand_biais(variant, th_sbiais, th_sb_cov, th_sb_pval):
            filters.append("StrandBiais")

        if vsup != float("nan") and rsup != float("nan") and (vsup / (vsup + rsup)) < rvt:
            filters.append("NoRealSupport")

        if len(filters) == 0:
            filters.append("PASS")

        variant.FILTER = filters

        writer.write_record(variant)


def compute_support(nodes, node2len, node_support):
    nodes = nodes.split("_")

    if len(nodes) <= 2:
        return float("nan")

    all_supports = 0
    path_len = 0

    for node in nodes[1:-1]:
        all_supports += max(node_support[node]["forward"]) + max(node_support[node]["reverse"])
        path_len += node2len[node]

    return all_supports / path_len


def strand_biais(record, th_sb_cov, th_sb_pval, th_sbiais):
    rcov = float(record.INFO["RCOV"])
    vcov = float(record.INFO["VCOV"])
    vcov_forward = float(record.INFO["VCOVF"])
    vcov_reverse = float(record.INFO["VCOVR"])
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
