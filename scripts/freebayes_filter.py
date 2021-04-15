import vcfpy

from shared import compute_entropy, alexSBFilter


def main(vcf, entropy_th, ratio_ref_cov_th, out):
    reader = vcfpy.Reader.from_path(vcf)

    writer = vcfpy.Writer.from_path(out, reader.header)

    for record in reader:
        coverage = record.INFO["DP"]
        ref_cov = record.INFO["RO"]
        ratio = ref_cov / coverage

        var_f = int(record.INFO["SAF"][0])
        var_r = int(record.INFO["SRF"])

        cov = var_f + var_r
        mincov = min(var_f, var_r)
        fq = mincov / cov
        minfq = min(1 - fq, fq)

        entropy = compute_entropy(record.REF)

        if (entropy < 0.001 or entropy > entropy_th) and ratio < ratio_ref_cov_th and not alexSBFilter(mincov, cov, minfq):

            writer.write_record(record)


if "snakemake" in locals():
    main(
        snakemake.input["vcf"],
        snakemake.params["e_th"],
        snakemake.params["ratio_ref_cov"],
        snakemake.output["output"],
    )
else:
    import sys

    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
