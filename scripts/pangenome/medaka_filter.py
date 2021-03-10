import vcfpy


def main(inputs, output):
    nano_var = None
    if "nanopolish" in inputs:
        reader = vcfpy.Reader.from_path(inputs["nanopolish"])
        for record in reader:
            key = (
                str(record.POS)
                + record.REF
                + "".join([a.serialize() for a in record.ALT])
            )
            record.INFO = dict()
            record.FORMAT = dict()
            nano_var[key] = record

    header = vcfpy.Header(samples=vcfpy.header.SamplesInfos(["all_var"]))
    # header = vcfpy.Header()
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    header.add_contig_line({"ID": "MN908947.3", "length": "29903"})
    header.add_filter_line({"Description": "All filters passed", "ID": "PASS"})
    header.add_filter_line({"Description": "dp", "ID": "dp"})

    writer = vcfpy.Writer.from_path(output, header)

    reader = vcfpy.Reader.from_path(inputs["medaka"])
    for vcf in reader:
        key = (
            str(record.POS)
            + record.REF
            + "".join([a.serialize() for a in record.ALT])
        )
        if any([alt.type == "DEL" for alt in record.ALT]) and nano_var is not None and key is not in nano_var:
            continue

        if record.INFO["AQ"] < config["medakaQualityThreshold"]:
            continue

        writer.write_record(record)


if "snakemake" in locals():
    main(snakemake.input, snakemake.output)
else:
    import sys

    main(sys.argv[1], sys.argv[2:])
