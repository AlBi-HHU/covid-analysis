
import vcfpy

from collections import defaultdict

def main(in_vcf, out_vcf):
    reader = vcfpy.Reader.from_path(in_vcf)

    header = reader.header
    
    pos2var = defaultdict(list)
    for record in reader:
        coverage = record.INFO["VCOV"] + record.INFO["RCOV"]
        pos2var[record.POS].append((coverage, record))

    header.add_info_line({"ID": "MULTIPLE", "Type": "Flag", "Number": "1", "Description": "Pangenome found multiple variant at this position"})
    writer = vcfpy.Writer.from_path(out_vcf, reader.header)
    
    for key, values in pos2var.items() :
        if len(values) == 1:
            record.INFO["MULTIPLE"] = False
            writer.write_record(values[0][1])
        else:
            record.INFO["MULTIPLE"] = True
            values.sort(key=lambda x: x[0], reverse=True)
            writer.write_record(values[0][1])


if "snakemake" in locals():
    main(snakemake.input[0], snakemake.output[0])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2])
