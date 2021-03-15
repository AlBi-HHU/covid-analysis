import vcfpy

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict

def main(in_vcf, min_cov, out_vcf):



    reader = vcfpy.Reader.from_path(in_vcf)

    header = reader.header
    
    pos2var = defaultdict(list)
    for record in reader:
        coverage = record.INFO["VCOV"] + record.INFO["RCOV"]
        pos2var[(record.POS, tuple(record.REF), tuple(record.ALT))].append((coverage, record))

    header.add_info_line({"ID": "MULTIPLE", "Type": "Flag", "Number": "1", "Description": "Pangenome found multiple variant at this position"})
    writer = vcfpy.Writer.from_path(out_vcf, reader.header)
    
    for key, values in pos2var.items():
        if len(values) == 1:
            coverage = values[0][0]
            variant = values[0][1]
        else:
            values.sort(key=lambda x: x[0], reverse=True)
            coverage = values[0][0]
            variant = values[0][1]
            variant.INFO["MULTIPLE"] = True

        if coverage > min_cov:
            writer.write_record(variant)


if "snakemake" in locals():
    main(snakemake.input['vcf'],snakemake.input['alignment'], snakemake.input['pangenome'],snakemake.params["min_cov"], snakemake.output[0])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4],sys.arg[5])
