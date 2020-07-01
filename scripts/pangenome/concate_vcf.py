
import vcfpy
import collections

def main(out_path, vcfs):

    header = vcfpy.Header(samples=vcfpy.header.SamplesInfos(["all_var"]))
    header.add_contig_line(vcfpy.OrderedDict({"ID": "MN908947.3", "length": "29903"}))

    writer = vcfpy.Writer.from_path(out_path, header)

    already_write = set()
    
    for vcf in vcfs:
        reader = vcfpy.Reader.from_path(vcf)

        for record in reader:
            h = str(record.POS) + record.REF + "".join([a.serialize() for a in record.ALT])
            if h in already_write:
                continue

            already_write.add(h)
            
            ## vcfpy black magic I didn't understand what I do but it's work don't touch it
            record.INFO = dict()
            record.FORMAT = dict()
            record.calls = [vcfpy.Call(sample="all_var", data=vcfpy.OrderedDict(), site=record)]
            record.call_for_sample["all_var"] = record.calls[0]
            ## end of black magic
            
            writer.write_record(record)

if "snakemake" in locals():
    main(snakemake.output["out"], snakemake.input["vcfs"])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2:])
