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

    reader = vcfpy.Reader.from_path(inputs["medaka"])

    writer = vcfpy.Writer.from_path(output, reader.header)

    for record in reader:
        key = (
            str(record.POS)
            + record.REF
            + "".join([a.serialize() for a in record.ALT])
        )

        if any([alt.type == "DEL" for alt in record.ALT]) and nano_var is not None and key not in nano_var:
            continue

        if record.INFO["AQ"] < snakemake.config["medakaQualityThreshold"]:
            continue

        writer.write_record(record)


if "snakemake" in locals():
    main(snakemake.input, snakemake.output)
else:
    import sys

    main(sys.argv[1], sys.argv[2:])
