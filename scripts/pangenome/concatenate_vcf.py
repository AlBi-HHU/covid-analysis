import os
import vcfpy


def main(out_path, vcfs):
    header = vcfpy.Header(samples=vcfpy.header.SamplesInfos(["all_var"]))
    # header = vcfpy.Header()
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    header.add_contig_line({"ID": "MN908947.3", "length": "29903"})
    header.add_filter_line({"Description": "All filters passed", "ID": "PASS"})
    header.add_filter_line({"Description": "dp", "ID": "dp"})
    header.add_info_line(
        {
            "ID": "ORI",
            "Number": 1,
            "Type": "String",
            "Description": "Where this variant come from if only first origin is peak if multiple origin exist",
        }
    )
    # header.add_line(vcfpy.SampleHeaderLine.from_mapping({'ID':'all_var'}))

    writer = vcfpy.Writer.from_path(out_path, header)

    # Keep track of where each variant comes from
    variant_sources = {}
    variants = {}

    # Pass 1: Identify Variants and Sources
    for vcf in vcfs:

        if os.stat(vcf).st_size == 0:  # filter for empty vcfs generated by ctx
            # TODO: Maybe log or notify
            continue

        try:
            reader = vcfpy.Reader.from_path(vcf)
            for record in reader:
                h = (
                    str(record.POS)
                    + "/"
                    + record.REF
                    + "".join([a.serialize() for a in record.ALT])
                )
                # If we see this the first time we save it and initialize the source list
                if h not in variants:
                    variant_sources[h] = set()
                    variants[h] = record
                variant_sources[h].add(get_origin(vcf))
        except:
            raise Exception("Error parsing vcf: {}: {}".format(vcf, sys.exc_info()[0]))

    # Pass 2: Write to output
    for variant in sorted(variants.keys(), key=lambda x: int(x.split("/")[0])):
        record = variants[variant]
        ## vcfpy black magic I didn't understand what I do but it's work don't touch it  // :)
        record.INFO = {"ORI": "/".join(list(variant_sources[variant]))}
        record.FORMAT = dict()

        record.calls = [
            vcfpy.Call(sample="all_var", data=vcfpy.OrderedDict(), site=record)
        ]
        record.call_for_sample["all_var"] = record.calls[0]

        ## end of black magic

        writer.write_record(record)


def get_origin(path):
    if "freebayes" in path:
        return "Freebayes"
    elif "arti" in path:  # Medaka is filtered one more
        return "Medaka"
    elif "nanopolish" in path:
        return "Nanopolish"
    else:
        return "GISAID"


if "snakemake" in locals():
    main(snakemake.output["out"], snakemake.input["vcfs"])
else:
    import sys

    main(sys.argv[1], sys.argv[2:])
