
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
        
        entropy = compute_entropy(record.REF)

        if (entropy < 0.001 or entropy > entropy_th) and ratio < ratio_ref_cov_th:
            
            writer.write_record(record)

            
def compute_entropy(seq):
    entropy = 0
    for k, v in Counter(seq).items():
        entropy += v/len(seq) * log2(v/len(seq))

    return entropy * -1

    
if "snakemake" in locals():
    main(snakemake.input["vcf"], snakemake.params["e_th"], snakemake.params["ratio_ref_cov"], snakemake.output["output"])
else:
    import sys
    
    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4])
