
import vcfpy

from math import log2
from collections import Counter


def main(vcf, entropy_th, ratio_ref_cov_th, out):
    reader = vcfpy.Reader.from_path(vcf)

    writer = vcfpy.Writer.from_path(out, reader.header)

    for record in reader:
        coverage = record.INFO["DP"]
        ref_cov = record.INFO["RO"]
        ratio = ref_cov / coverage

        seqs = [record.REF]
        seqs.extend([alt.serialize() for alt in record.ALT])
        homo = any([is_homopolymer(seq) for seq in seqs])

        if not homo and ratio < ratio_ref_cov_th:
            writer.write_record(record)


def is_homopolymer(seq):
    successive_base = 1
    for b, c in zip(seq, seq[1:]):
        if b == c:
            successive_base += 1
        else:
            successive_base = 1

        if successive_base >= 4:
            return True

    return False


if "snakemake" in locals():
    main(snakemake.input["vcf"], snakemake.params["e_th"], snakemake.params["ratio_ref_cov"], snakemake.output["output"])
else:
    import sys

    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
