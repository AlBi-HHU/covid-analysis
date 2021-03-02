import os
import vcfpy

def main(out_path, vcfs):

    header = vcfpy.Header(samples=vcfpy.header.SamplesInfos(["all_var"]))
    #header = vcfpy.Header()
    header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv4.2'))
    header.add_contig_line({"ID": "MN908947.3", "length": "29903"})
    header.add_filter_line({'Description': "All filters passed", 'ID': 'PASS'})
    header.add_filter_line({'Description': "dp", 'ID': 'dp'})
    #header.add_line(vcfpy.SampleHeaderLine.from_mapping({'ID':'all_var'}))
    
    writer = vcfpy.Writer.from_path(out_path, header)

    already_write = set()
    
    for vcf in vcfs:

        if os.stat(vcf).st_size == 0: #filter for empty vcfs generated by ctx
            #TODO: Maybe log or notify
            continue

        try:
            reader = vcfpy.Reader.from_path(vcf)
            for record in reader:
                h = str(record.POS) + record.REF + "".join([a.serialize() for a in record.ALT])
                if h in already_write:
                    continue

                already_write.add(h)
                ## vcfpy black magic I didn't understand what I do but it's work don't touch it  // :)
                record.INFO = dict()
                record.FORMAT = dict()


                record.calls = [vcfpy.Call(sample="all_var", data=vcfpy.OrderedDict(), site=record)]
                record.call_for_sample["all_var"] = record.calls[0]


                ## end of black magic

                writer.write_record(record)
        except:
            raise Exception('Error parsing vcf: {}: {}'.format(vcf,sys.exc_info()[0]))
            


if "snakemake" in locals():
    main(snakemake.output["out"], snakemake.input["vcfs"])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2:])
