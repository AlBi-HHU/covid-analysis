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

    '''
    if record.FILTER != ["PASS"]:
        continue
    '''

    # We only have single variants
    ref = record.REF
    alt = record.ALT[
        0
    ].value  # therefore picking the first one picks the only existing variant
    pos = record.POS

    if pos in pileup:
        logging.debug("Corresponding pileup record: {}".format(pileup[pos]))

    if 'LowRVT' in record.FILTER:
        continue

    #Check all possible filter flags
    RCOVERAGE = 'RCoverage' in record.FILTER
    VCOVERAGE = 'VCoverage' in record.FILTER
    RSTRANDBIAS = 'RStrandBias' in record.FILTER
    VSTRANDBIAS = 'VStrandBias' in record.FILTER

    mask = str(int(RCOVERAGE))+str(int(RSTRANDBIAS))+str(int(VCOVERAGE))+str(int(VSTRANDBIAS))

    if mask == '1111':
        raise Exception('Dropout, should have been caught before!')
    elif mask == '1011':
        raise Exception('Dropout, should have been caught before!')
    elif mask == '1110':
        raise Exception('Dropout, should have been caught before!')
    elif mask == '1010':
        raise Exception('Dropout, should have been caught before!')
    elif mask == '0000':
        varRatio = float(record.INFO["CORHETRATIO"])

        if th_het <= varRatio <= 1 - th_het:
            if len(alt) == 1 and len(ref) == 1:

                record.ALT[0].value = ambiguityLetters[
                    frozenset({record.ALT[0].value, record.REF})
                ]
            else:
                # Here, we have a heterozygous structural variant, since we can't really encode this in a linear consensus we just keep the information in the .vcf file
                record.INFO["HSV"] = True
    elif mask == '0101':
        varRatio = sum(record.INFO['VCOV'])/sum(record.INFO['RCOV'])
        if th_het <= varRatio <= 1 - th_het:
            if len(alt) == 1 and len(ref) == 1:

                record.ALT[0].value = ambiguityLetters[
                    frozenset({record.ALT[0].value, record.REF})
                ]
            else:
                # Here, we have a heterozygous structural variant, since we can't really encode this in a linear consensus we just keep the information in the .vcf file
                record.INFO["HSV"] = True
    elif mask == '0001':
        continue #REF
    elif mask == '0010':
        continue #REF
    elif mask == '0011':
        continue #REF
    elif mask == '0110':
        continue #REF
    elif mask == '0111':
        continue #REF
    elif mask == '1000':
        pass #HOM ALT
    elif mask == '1001':
        pass #HOM ALT
    elif mask == '1100':
        pass #HOM ALT
    elif mask == '1101':
        pass #HOM ALT
    elif mask == '0100':
        pass #HOM ALT
    record.FILTER = ['PASS']
    writer.write_record(record)
