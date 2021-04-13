import vcfpy

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure

from shared import (
    get_vcf_strand_bias_filter,
    parsePileupStrandAwareLight,
    ambiguityLetters,
)
import logging

# Enable logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)

# Reassign parameters from snakemake config for better readability and type casting (?)
th_cov = int(snakemake.config["consensusMinCov"])
th_het = float(snakemake.config["thresholdHomCall"])

# See shared.py
pileup = parsePileupStrandAwareLight(snakemake.input["pileup"])

reader = vcfpy.Reader.from_path(snakemake.input["vcf"])

# add additional header line to mark het vars
header = reader.header  # Copy the existing header
header.add_info_line(
    {
        "ID": "HSV",
        "Type": "Flag",
        "Number": "1",
        "Description": "Variant might be a heterozygous SV",
    }
)


# First of all we identify drop-out regions where our coverage is too low to make any calls, we will mask them with N letters
with open(snakemake.output["nMask"], "w") as outfile:
    for pos in range(
        1, snakemake.config["ref_genome_length"] + 1
    ):  # Note, that coordinates are 1-based
        # Positions that are not in the pileup are not covered at all, for all others: Check the coverage and compare to the threshold
        if (not (pos in pileup)) or (sum(pileup[pos].values()) < th_cov):
            outfile.write("{}\t{}\n".format(snakemake.config["ref_genome_chr"], pos))

writer = vcfpy.Writer.from_path(snakemake.output["vcf"], header)

for record in reader:
    logging.debug("Processing record: {}".format(record))

    if record.FILTER != ["PASS"]:
        continue

    # We only have single variants
    ref = record.REF
    alt = record.ALT[
        0
    ].value  # therefore picking the first one picks the only existing variant
    pos = record.POS

    if pos in pileup:
        logging.debug("Corresponding pileup record: {}".format(pileup[pos]))

    varRatio = float(record.INFO["CORHETRATIO"])

    if th_het <= varRatio <= 1 - th_het:

        sb_ref = get_vcf_strand_bias_filter(
            record, "R"
        )  # True if there is a strand bias issue in the reference component of the HET call -> Assume a pure HOM call in this case
        if sb_ref:
            # Nothing needs to be done here: We will simply not substitute the matching ambiguity character, thus "keeping" the pure HOM call
            pass
        # SNPs get the ambiguous base characterss
        elif len(alt) == 1 and len(ref) == 1:

            record.ALT[0].value = ambiguityLetters[
                frozenset({record.ALT[0].value, record.REF})
            ]
        else:
            # Here, we have a heterozygous structural variant, since we can't really encode this in a linear consensus we just keep the information in the .vcf file
            record.INFO["HSV"] = True

    writer.write_record(record)
