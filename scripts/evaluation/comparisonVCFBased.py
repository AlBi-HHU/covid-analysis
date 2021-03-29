import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import (
    parsePileupStrandAwareLight,
    getTotalCoverage,
    ambiguityLetters,
    ambiguityLetters_inverted,
    getCoverage,
    getMinorStrandAbs,
    getMinorStrandFrequency,
    alexSBFilter,
    isAmbiguous,
)
from Bio import SeqIO
import vcfpy
import pandas as pd

### Skript to perform a vcf based comparison between variant calls


# Counter Variables for Top Level Stats
cnt_realVariants = 0  #

cnt_realHETSNPs = 0  #
cnt_detectedHETSNPs = 0  #

cnt_realSNP = 0  #
cnt_detectedSNP = 0  #

cnt_realINS = 0  #
cnt_detectedINS = 0  #

cnt_realHETINS = 0
cnt_detectedHETINS = 0

cnt_realDEL = 0  #
cnt_detectedDEL = 0  #

cnt_realHETDEL = 0
cnt_detectedHETDEL = 0

cnt_concordance = 0  #
cnt_falsePositives = 0
cnt_falseNegatives = 0
cnt_discordance = 0  #
cnt_detectedVariants = 0  #
cnt_relevantPositions = 0  #
cnt_comparablePositions = 0  #
cnt_unscoredPositions = 0  #
cnt_illuminaDropouts = 0  #
cnt_nanoporeDropouts = 0  #

cnt_illuminaDropouts_unscored = 0  #
cnt_nanoporeDropouts_unscored = 0  #

cnt_implicit_agreement = 0  #

cnt_unfairComparisons = 0  #

cnt_falseHomozygous = 0  #

# TODO: Track structural SV, Track detected Heterozygosity based on Medaka / Nano Info

# We keep track of tuples in addition to construct a pandas dataframe for easy transformation into altair plots
dataTuples = []

# Keep Output files open
with open(snakemake.output["text"], "w") as outfile, open(
    snakemake.output["filter"], "w"
) as filterfile, open(snakemake.log[0], "w") as logfile:

    # Read required input
    reference = SeqIO.read(snakemake.input["ref"], "fasta")

    # Read Illumina Average Coverage Track
    df = pd.read_csv(snakemake.input["cov"])
    averageIlluminaCoverage = {}
    for row in df.itertuples():
        averageIlluminaCoverage[row.pos] = row.cov

    # Input comes in blocks of fours
    data = iter(snakemake.input["comparisonFiles"])
    for ivarPseudoVCF in data:

        fid = ivarPseudoVCF.split("/")[-2]

        # Write down file name to keep output file readable
        outfile.write(fid + "\n")

        # Variable Reassignment to keep code more readable
        ivarPseudoVCF = pd.read_csv(ivarPseudoVCF, sep="\t")
        # ivarPseudoVCF = ivarPseudoVCF[ivarPseudoVCF.PASS != False]  # Filter passed vars only

        pancovVCF = vcfpy.Reader.from_path(next(data))

        illuminapileup = parsePileupStrandAwareLight(next(data))
        nanoporepileup = parsePileupStrandAwareLight(next(data))

        ### Step 1: Determine all relevant positions
        relevantPositions = set()

        recordsNanopore = {}
        recordsIllumina = {}
        logfile.write("processing {}".format(fid))
        for record in pancovVCF:
            relevantPositions.add(record.POS)
            recordsNanopore[record.POS] = record

        for position in ivarPseudoVCF["POS"].unique():

            record = ivarPseudoVCF[ivarPseudoVCF.POS == position]
            altallele = record["ALT"].values[0]
            # if we have a deletion or an insertion we ignore it for the sb test
            if altallele.startswith("-") or altallele.startswith("+"):
                recordsIllumina[position] = record
                relevantPositions.add(position)
                continue  # with the next position

            # otherwise we apply the strand bias filter test
            components = (
                ambiguityLetters_inverted[altallele]
                if altallele in ambiguityLetters_inverted
                else [altallele]
            )

            for component in components:

                cov = getCoverage(illuminapileup[position], altallele)
                abs = getMinorStrandAbs(illuminapileup[position], altallele)
                fq = getMinorStrandFrequency(illuminapileup[position], altallele)

                if alexSBFilter(cov, abs, fq):
                    filterfile.write(
                        "{}: iVar call for {} did not pass Alex Filter with AltAllele Cov {}: Total Cov: {} and Minor Strand FQ: {}\n".format(
                            fid, altallele, cov, abs, fq
                        )
                    )
                    break
                elif str(record["PASS"].values[0]) == "FALSE":
                    filterfile.write(
                        "{}: iVar call for {} did not pass the internal iVar filter\n".format(
                            fid, altallele
                        )
                    )
                    break
            else:
                recordsIllumina[position] = record
                relevantPositions.add(position)

        ### Step 2: Process

        fields = ["pos", "ref", "ivar", "nanopore-method"]
        outfile.write("\t".join(fields) + "\n")

        # Loop over the entire genome
        for position in range(1, snakemake.config["ref_genome_length"] + 1):
            # Determine Nanopore and Illumina Coverage
            nanoporeCoverage = (
                getTotalCoverage(nanoporepileup[position])
                if position in nanoporepileup
                else 0
            )
            illuminaCoverage = (
                getTotalCoverage(illuminapileup[position])
                if position in illuminapileup
                else 0
            )

            # Decide whether the position counts or not
            nanoporeDropout = (
                nanoporeCoverage < snakemake.config["nanoporeCoverageCutoff"]
            )

            illuminaDropout = illuminaCoverage < (
                averageIlluminaCoverage[position]
                * snakemake.config["illuminaCoverageCutoff"]
            )

            # Check non-relevant positions
            if position not in relevantPositions:
                cnt_illuminaDropouts_unscored += illuminaDropout
                cnt_nanoporeDropouts_unscored += nanoporeDropout
                cnt_implicit_agreement += (not illuminaDropout) and (
                    not nanoporeDropout
                )
                continue  # with next position
            else:
                # From here on: relevant position

                # Track some properties for easier visualization of results later on
                bool_falseNegative = False
                bool_falsePositive = False
                bool_heterozygousIllu = False
                bool_heterozygousNano = False
                bool_concordance = False

                # Dropouts
                if illuminaDropout or nanoporeDropout:
                    cnt_unscoredPositions += 1
                    if illuminaDropout:
                        cnt_illuminaDropouts += 1
                    if nanoporeDropout:
                        cnt_nanoporeDropouts += 1

                # Resolve Type and Value of Variant

                nanoporeType = "Not Called"  # INS,DEL,SNP
                nanoporeValue = ""  # Length of Del / Alt Allele / Insertion Seq

                if position in recordsNanopore:
                    nanoporeType = recordsNanopore[position].ALT[0].type
                    if nanoporeType == "INS":
                        nanoporeValue = (
                            recordsNanopore[position].ALT[0].value[2:]
                        )  # Ignore first char as this is REF
                        cnt_detectedINS += 1
                        if (
                            snakemake.params["method"] == "pancov"
                            and "HSV" in recordsNanopore[position].INFO
                        ):
                            bool_heterozygousNano = True
                            cnt_detectedHETINS += 1
                            nanoporeValue += "(HET)"
                    elif nanoporeType == "DEL":
                        nanoporeValue = str(
                            len(recordsNanopore[position].REF) - 1
                        )  # Ignore first char as this is retained
                        cnt_detectedDEL += 1
                        if (
                            snakemake.params["method"] == "pancov"
                            and "HSV" in recordsNanopore[position].INFO
                        ):
                            bool_heterozygousNano = True
                            cnt_detectedHETDEL += 1
                            nanoporeValue += "(HET)"
                    elif nanoporeType == "SNV":
                        nanoporeValue = recordsNanopore[position].ALT[0].value
                        cnt_detectedSNP += 1
                        if snakemake.params["method"] == "pancov" and isAmbiguous(
                            nanoporeValue
                        ):
                            bool_heterozygousNano = True
                            cnt_detectedHETSNPs += 1

                illuminaType = "Not Called"  # INS,DEL,SNP
                illuminaValue = ""  # Length of Del / Alt Allele / Insertion Seq

                if position in recordsIllumina:
                    altval = recordsIllumina[position]["ALT"].values[0]
                    refval = recordsIllumina[position]["REF"].values[0]
                    altfreq = float(recordsIllumina[position]["ALT_FREQ"].values[0])
                    if altval.startswith("-"):
                        illuminaType = "DEL"
                        cnt_realDEL += 1
                    elif altval.startswith("+"):
                        illuminaType = "INS"
                        cnt_realINS += 1
                    else:
                        illuminaType = "SNV"
                        cnt_realSNP += 1

                    if illuminaType == "INS":
                        illuminaValue = altval[2:]  # Don't use the +
                        if (
                            snakemake.config["thresholdHomCall"]
                            <= altfreq
                            <= (1 - snakemake.config["thresholdHomCall"])
                        ):
                            bool_heterozygousIllu = True
                            cnt_realHETINS += 1
                            illuminaValue += (
                                "(HET)"  # Add (HET) as marker for heterozygosity
                            )
                    elif illuminaType == "DEL":
                        illuminaValue = str(
                            len(altval) - 1
                        )  # Ignore first char as this is retained
                        if (
                            snakemake.config["thresholdHomCall"]
                            <= altfreq
                            <= (1 - snakemake.config["thresholdHomCall"])
                        ):
                            bool_heterozygousIllu = True
                            cnt_realHETDEL += 1
                            illuminaValue += "(HET)"
                    elif illuminaType == "SNV":
                        illuminaValue = altval
                        if (
                            snakemake.config["thresholdHomCall"]
                            <= altfreq
                            <= (1 - snakemake.config["thresholdHomCall"])
                        ):
                            bool_heterozygousIllu = True
                            cnt_realHETSNPs += 1
                            illuminaValue = ambiguityLetters[
                                frozenset((altval, refval))
                            ]

                if not (illuminaDropout or nanoporeDropout):
                    if (illuminaType == nanoporeType) and (
                        illuminaValue == nanoporeValue
                    ):
                        bool_concordance = True
                        cnt_concordance += 1
                    else:
                        cnt_discordance += 1
                        if (position in recordsIllumina) and (
                            position not in recordsNanopore
                        ):
                            bool_falseNegative = True
                            cnt_falseNegatives += 1
                        if (position in recordsNanopore) and (
                            position not in recordsIllumina
                        ):
                            bool_falsePositive = True
                            cnt_falsePositives += 1
                        if bool_heterozygousIllu and not bool_heterozygousNano:
                            cnt_falseHomozygous += 1
                        if (
                            bool_heterozygousIllu
                            and snakemake.params["method"] != "pancov"
                        ):
                            cnt_unfairComparisons += 1

            # Write Text Output
            outfile.write(
                "{}\t{}\t{}\t{}\n".format(
                    position,
                    reference[int(position - 1)],  # SeqIO is 0-based
                    illuminaType + " " + illuminaValue,
                    nanoporeType + " " + nanoporeValue,
                )
            )

            # Write to tuples
            dataTuples.append(
                (
                    fid,
                    position,
                    reference[position - 1],
                    illuminaType,
                    illuminaValue,
                    nanoporeType,
                    nanoporeValue,
                    nanoporeDropout,
                    illuminaDropout,
                    bool_falseNegative,
                    bool_falsePositive,
                    bool_heterozygousIllu,
                    bool_heterozygousNano,
                    bool_concordance,
                    list(illuminapileup[position] for position in sorted(illuminapileup))
                    if position in illuminapileup
                    else "Dropout",
                    list(nanoporepileup[position] for position in sorted(nanoporepileup))
                    if position in nanoporepileup
                    else "Dropout",
                    reference[position - 1-3:position-1+3+1].seq, #3 preceding, 3 succeeding bases
                )
            )

        # Calculate some additional stats
        cnt_detectedVariants += len(recordsNanopore)
        cnt_realVariants += len(recordsIllumina)
        cnt_relevantPositions += len(relevantPositions)

    cnt_comparablePositions = cnt_relevantPositions - cnt_unscoredPositions

    outfile.write("Real (iVar) Variants:{} \n".format(cnt_realVariants))

    outfile.write("Real (iVar) SNPs:{} \n".format(cnt_realSNP))
    outfile.write("Detected SNPs:{} \n".format(cnt_detectedSNP))
    outfile.write("Real (iVar) HET SNPs:{} \n".format(cnt_realHETSNPs))
    outfile.write("Detected HET SNPs:{} \n".format(cnt_detectedHETSNPs))

    outfile.write("Real (iVar) INS:{} \n".format(cnt_realINS))
    outfile.write("Detected INS:{} \n".format(cnt_detectedINS))
    outfile.write("Real (iVar) HET INS:{} \n".format(cnt_realHETINS))
    outfile.write("Detected HET INS:{} \n".format(cnt_detectedHETINS))

    outfile.write("Real (iVar) DEL:{} \n".format(cnt_realDEL))
    outfile.write("Detected DEL:{} \n".format(cnt_detectedDEL))
    outfile.write("Real (iVar) HET DEL:{} \n".format(cnt_realHETDEL))
    outfile.write("Detected HET DEL:{} \n".format(cnt_detectedHETDEL))

    outfile.write(
        "Relevant Positions: {} (of which {} could not be evaluated)\n".format(
            cnt_relevantPositions, cnt_unscoredPositions
        )
    )
    outfile.write(
        "Concordance:{} of {} comparable positions \n".format(
            cnt_concordance, cnt_comparablePositions
        )
    )
    outfile.write("FP:{} \n".format(cnt_falsePositives))
    outfile.write("FN:{} \n".format(cnt_falseNegatives))
    outfile.write(
        "Discordance: {} of {} comparable positions \n".format(
            cnt_discordance, cnt_comparablePositions
        )
    )
    outfile.write(
        "Possible unfair comparisons for method != pancov (HET was not found or only as HOM): {} \n".format(
            cnt_unfairComparisons
        )
    )
    outfile.write(
        "Positions where a HOM variant was detected but iVar detected a HET: {} \n".format(
            cnt_falseHomozygous
        )
    )
    outfile.write("Detected Variants: {} \n".format(cnt_detectedVariants))
    outfile.write(
        "Unscored Positions: {} ({} Illumina and {} Nanopore Dropouts) \n".format(
            cnt_unscoredPositions, cnt_illuminaDropouts, cnt_nanoporeDropouts
        )
    )

    # Calculate top level stats
    precision = cnt_concordance / (cnt_concordance + cnt_falsePositives)
    recall = cnt_concordance / (cnt_concordance + cnt_falseNegatives)
    f1 = (2 * cnt_concordance) / (
        2 * cnt_concordance + cnt_falsePositives + cnt_falseNegatives
    )
    outfile.write("Precision (Concordance / (Concordance+FP)): {}\n".format(precision))
    outfile.write("Recall (Concordance / (Concordance+FN)): {}\n".format(recall))
    outfile.write("F1 ( (2*Concordance)/(2*Concordance+FP+FN)): {}\n".format(f1))

    # Additional TLS
    accuracy = cnt_concordance / cnt_comparablePositions
    outfile.write("Accuracy (Concordance/Comparable Positions): {}\n".format(accuracy))
    outfile.write(
        "Additional Nanopore Dropouts without VCs in either method: {}\n".format(
            cnt_nanoporeDropouts_unscored
        )
    )
    outfile.write(
        "Additional Illumina Dropouts without VCs in either method: {}\n".format(
            cnt_illuminaDropouts_unscored
        )
    )
    outfile.write(
        "Implicit reference agreement (IA): {}\n".format(cnt_implicit_agreement)
    )

    # Calculate top level stats in a variation
    precision_var = (cnt_concordance + cnt_implicit_agreement) / (
        cnt_concordance + cnt_falsePositives + cnt_implicit_agreement
    )
    outfile.write(
        "Precision (Concordance+IA / (Concordance+FP+IA)): {}\n".format(precision)
    )
    accuracy = (cnt_concordance + cnt_implicit_agreement) / (
        cnt_comparablePositions + cnt_implicit_agreement
    )
    outfile.write(
        "Accuracy ((Concordance+IA)/(Comparable Positions+IA)): {}\n".format(accuracy)
    )


# Write Pandas Dataframe
df = pd.DataFrame(
    dataTuples,
    columns=[
        "sample",
        "position",
        "reference",
        "illuminaType",
        "illuminaValue",
        "nanoporeType",
        "nanoporeValue",
        "nanoporeDropout",
        "illuminaDropout",
        "falseNegative",
        "falsePositive",
        "realHET",
        "detectedHET",
        "concordance",
        "illuminaPileup",
        "nanoporePileup",
        "referenceContext"
    ],
)

df.to_csv(snakemake.output["dataframe"])
