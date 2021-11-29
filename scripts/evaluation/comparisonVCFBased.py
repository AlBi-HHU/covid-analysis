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

#cnt_concordance_HOM = 0  #
#cnt_concordance_HET = 0
#cnt_falsePositives_HOM = 0
#cnt_falsePositives_HET = 0
#cnt_falseNegatives_HOM = 0
#cnt_falseNegatives_HET = 0
#cnt_discordance = 0  #

cnt_detectedVariants = 0  #
cnt_relevantPositions = 0  #
cnt_comparablePositions = 0  #
cnt_unscoredPositions = 0  #
cnt_illuminaDropouts = 0  #
cnt_nanoporeDropouts = 0  #


cnt_tp_perallele_het = 0
cnt_fp_perallele_het = 0
cnt_fn_perallele_het = 0

cnt_tp_perallele_hom = 0
cnt_fp_perallele_hom = 0
cnt_fn_perallele_hom = 0

cnt_con_perposition_hom = 0
cnt_dis_perposition_hom = 0

cnt_con_perposition_het = 0
cnt_dis_perposition_het = 0

cnt_illuminaDropouts_unscored = 0  #
cnt_nanoporeDropouts_unscored = 0  #

cnt_implicit_agreement = 0  #

cnt_unfairComparisons = 0  #

cnt_falseHomozygous = 0  #

cnt_multiallelic_tp = 0
cnt_multiallelic_fn = 0
cnt_multiallelic_fp = 0
cnt_multiallelic_tn = 0

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

            # We keep track of components that are either filtered or confirmed as true
            # For example if we have Y we keep track of C and T separately and we might remove one component
            filteredComponents = []

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
                    filteredComponents.append(component)
                elif str(record["PASS"].values[0]) == "FALSE":
                    filterfile.write(
                        "{}: iVar call for {} did not pass the internal iVar filter\n".format(
                            fid, altallele
                        )
                    )
                    filteredComponents.append(component)

            # Remove the filtered Components
            components = set(components) - set(filteredComponents)

            # Check for two cases:
            # Case 1: The only remaining component is REF
            if len(components) == 1 and list(components)[0] == record["REF"].values[0]:
                pass
            # Case 2: All components got filtered by SB
            elif len(components) == 0:
                pass
            else:
                # Determine the new alt value depending on the filtering done before
                record["ALT"].values[0] = ambiguityLetters[frozenset(components)]
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
                        elif snakemake.params["method"] == "nanopolish" and (
                            snakemake.config["thresholdHomCall"]
                            <= recordsNanopore[position].INFO["SupportFraction"]
                            <= (1 - snakemake.config["thresholdHomCall"])
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
                        elif snakemake.params["method"] == "nanopolish" and (
                            snakemake.config["thresholdHomCall"]
                            <= recordsNanopore[position].INFO["SupportFraction"]
                            <= (1 - snakemake.config["thresholdHomCall"])
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
                        elif snakemake.params["method"] == "nanopolish" and (
                            snakemake.config["thresholdHomCall"]
                            <= recordsNanopore[position].INFO["SupportFraction"]
                            <= (1 - snakemake.config["thresholdHomCall"])
                        ):
                            bool_heterozygousNano = True
                            cnt_detectedHETSNPs += 1
                            altval = recordsNanopore[position].ALT[0].value
                            refval = recordsNanopore[position].REF
                            nanoporeValue = ambiguityLetters[
                                frozenset((altval, refval))
                            ]
                    else:
                        print('Case not handled: {} / {}'.format(nanoporeType,nanoporeValue))
                        assert False
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

                ################
                #S C O R I N G
                #
                ################

                if not (illuminaDropout or nanoporeDropout):

                    #Multi-Allelic Classification
                    if bool_heterozygousNano and not bool_heterozygousIllu:
                        cnt_multiallelic_fp += 1
                    elif bool_heterozygousIllu and not bool_heterozygousNano:
                        cnt_multiallelic_fn += 1
                    elif bool_heterozygousIllu and bool_heterozygousNano:
                        cnt_multiallelic_tp += 1
                    elif not bool_heterozygousIllu and not bool_heterozygousNano:
                        cnt_multiallelic_tn += 1

                    if (
                        (illuminaType == nanoporeType)
                        and (illuminaValue == nanoporeValue)
                        and (bool_heterozygousIllu == bool_heterozygousNano)
                    ):
                        # if this is a het position
                        if bool_heterozygousNano:
                            # Both, alt and ref got called correctly
                            cnt_tp_perallele_het += 2
                            cnt_con_perposition_het += 1
                        else:
                            cnt_tp_perallele_hom += 1
                            cnt_con_perposition_hom += 1
                        bool_concordance = True
                    else: #Discordance, there is some sort of disagreement
                        #cnt_discordance += 1
                        if ( #clean, missed HOM call
                            (position in recordsIllumina)
                            and (position not in recordsNanopore)
                            and (not bool_heterozygousIllu)
                        ):
                            bool_falseNegative = True
                            cnt_fn_perallele_hom += 1
                            cnt_dis_perposition_hom += 1
                            #cnt_falseNegatives_HOM += 1
                        if ( #clean, false positive HOM call
                            (position in recordsNanopore)
                            and (position not in recordsIllumina)
                            and (not bool_heterozygousNano)
                        ):
                            bool_falsePositive = True
                            #cnt_falsePositives_HOM += 1
                            cnt_fp_perallele_hom += 1
                            cnt_dis_perposition_hom += 1

                        if bool_heterozygousIllu and not bool_heterozygousNano:
                            cnt_falseHomozygous += 1

                        if (
                            bool_heterozygousIllu
                        ):
                            cnt_unfairComparisons += 1

                        # HET Stats
                        # Case 1: It's a real HET but we don't find it at all
                        if (
                            (position in recordsIllumina)
                            and (position not in recordsNanopore)
                            and bool_heterozygousIllu
                        ):
                            cnt_tp_perallele_het += 1
                            cnt_fn_perallele_het += 1
                            cnt_dis_perposition_het += 1
                            #cnt_falseNegatives_HET += 1
                            #cnt_concordance_HET += 1
                            bool_falseNegative = True
                        # Case 2: It's a called HET but Illumina knows of no variant at this location
                        if (
                            (position in recordsNanopore)
                            and (position not in recordsIllumina)
                            and bool_heterozygousNano
                        ):
                            cnt_tp_perallele_het += 1 #called the ref part correctly
                            cnt_fp_perallele_het += 1 #called the variant part incorrectly
                            cnt_dis_perposition_het += 1
                            #cnt_falsePositives_HET += 1
                            bool_falsePositive = True
                            #cnt_concordance_HET += 1
                        # Case 3: It's called in both cases but there is disagreement
                        # Case 3.1: het call in illu and hom call in nano
                        if (bool_heterozygousIllu) and (position in recordsNanopore):
                            if illuminaType == "SNV":
                                cnt_tp_perallele_het += 1 #variant part called correctly
                                cnt_fn_perallele_het += 1 #ref part called incorrectly
                                cnt_dis_perposition_het += 1
                                #cnt_falseNegatives_HET += 1
                                #cnt_concordance_HET += 1
                                bool_falseNegative = True
                            else:
                                if illuminaValue == nanoporeValue:
                                    '''
                                    bool_falseNegative = True
                                    cnt_falseNegatives_HET += 1
                                    cnt_concordance_HET += 1
                                    '''
                                    cnt_tp_perallele_het += 1 #variant part called correctly
                                    cnt_fn_perallele_het += 1 #ref part called incorrectly
                                    cnt_dis_perposition_het += 1
                                    #cnt_falseNegatives_HET += 1
                                    #cnt_concordance_HET += 1
                                    bool_falseNegative = True
                        # Case 3.2 het call in nano and hom call in illu
                        if (bool_heterozygousNano) and (position in recordsIllumina):
                            if illuminaType == "SNV":
                                cnt_fp_perallele_het += 1
                                cnt_tp_perallele_het += 1
                                cnt_dis_perposition_het += 1
                                #bool_falsePositive = True
                                #cnt_falsePositives_HET += 1
                                #cnt_concordance_HET += 1
                            else:
                                if illuminaValue == nanoporeValue:
                                    cnt_fp_perallele_het += 1
                                    cnt_tp_perallele_het += 1
                                    cnt_dis_perposition_het += 1
                                    #bool_falsePositive = True
                                    #cnt_falsePositives_HET += 1
                                    #cnt_concordance_HET += 1

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
                    list(
                        "{}:{}".format(allele, illuminapileup[position][allele])
                        for allele in sorted(illuminapileup[position])
                    )
                    if position in illuminapileup
                    else "Dropout",
                    list(
                        "{}:{}".format(allele, nanoporepileup[position][allele])
                        for allele in sorted(nanoporepileup[position])
                    )
                    if position in nanoporepileup
                    else "Dropout",
                    reference[
                        position - 1 - 3 : position - 1 + 3 + 1
                    ].seq,  # 3 preceding, 3 succeeding bases
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
        "Detected Multi Allelics total : {} \n".format(
            cnt_detectedHETSNPs + cnt_detectedHETINS + cnt_detectedHETDEL
        )
    )

    outfile.write("Relevant Positions: {}\n".format(cnt_relevantPositions))
    outfile.write("Unscored Positions: {}\n".format(cnt_unscoredPositions))
    outfile.write("Illumina Dropouts: {}\n".format(cnt_illuminaDropouts))
    outfile.write("Nanopore Dropouts: {}\n".format( cnt_nanoporeDropouts))

    outfile.write("Per Position Stats: \n")

    #Aggregate sums
    cnt_con_perposition = cnt_con_perposition_hom + cnt_con_perposition_het
    cnt_dis_perposition = cnt_dis_perposition_hom + cnt_dis_perposition_het

    outfile.write('Concordance (full agreement) per Position: {}\n'.format(cnt_con_perposition))
    outfile.write('Concordance for homozygous variants : {}\n'.format(cnt_con_perposition_hom))
    outfile.write('Concordance for heterozygous variants: {}\n'.format(cnt_con_perposition_het))

    outfile.write('Discordance (disagreement in at least one allele) per Position: {}\n'.format(cnt_dis_perposition))
    outfile.write('Discordance for homozygous variants (no method called a HET variant) : {}\n'.format(cnt_dis_perposition_hom))
    outfile.write('Discordance for heterozygous variants (at least one method called a HET variant): {}\n'.format(cnt_dis_perposition_het))

    #outfile.write("HOM Concordance:{} \n".format(cnt_concordance_HOM))
    #outfile.write("HOM FP:{} \n".format(cnt_falsePositives_HOM))
    #outfile.write("HOM FN:{} \n".format(cnt_falseNegatives_HOM))

    #Aggregate Sums
    cnt_tp_perallele = cnt_tp_perallele_hom + cnt_tp_perallele_het
    cnt_fp_perallele = cnt_fp_perallele_hom +cnt_fp_perallele_het
    cnt_fn_perallele = cnt_fn_perallele_hom +cnt_fn_perallele_het

    outfile.write("Per Allele Stats: \n")

    outfile.write('True Positives: {}\n'.format(cnt_tp_perallele))
    outfile.write('True Positives, where no method made a multi-allelic call: {}\n'.format(cnt_tp_perallele_hom))
    outfile.write('True Positives, where at least one method made a multi-allelic call: {}\n'.format(cnt_tp_perallele_het))

    outfile.write('False Positives: {}\n'.format(cnt_fp_perallele))
    outfile.write('False Positives, where no method made a multi-allelic call: {}\n'.format(cnt_fp_perallele_hom))
    outfile.write('False Positives, where at least one method made a multi-allelic call: {}\n'.format(cnt_fp_perallele_het))

    outfile.write('False Negatives: {}\n'.format(cnt_fn_perallele))
    outfile.write('False Negatives, where no method made a multi-allelic call: {}\n'.format(cnt_fn_perallele_hom))
    outfile.write('False Negatives, where at least one method made a multi-allelic call: {}\n'.format(cnt_fn_perallele_het))

    outfile.write('True Negatives are of infinite size due to an theoretically arbitrary amount of alternative alleles \n')

    outfile.write(
        "Possible unfair comparisons for methods that don't detect multi-allelic (multi-allelic was not found or only as HOM): {} \n".format(
            cnt_unfairComparisons
        )
    )
    outfile.write(
        "Positions where a HOM variant was detected but iVar detected a multi-allelic: {} \n".format(
            cnt_falseHomozygous
        )
    )
    outfile.write("Detected Variants: {} \n".format(cnt_detectedVariants))


    precision_HOM = cnt_tp_perallele_hom / (cnt_tp_perallele_hom + cnt_fp_perallele_hom)
    recall_HOM = cnt_tp_perallele_hom / (cnt_tp_perallele_hom + cnt_fn_perallele_hom)
    f1_HOM = (2*cnt_tp_perallele_hom) / (2 * cnt_tp_perallele_hom + cnt_fp_perallele_hom + cnt_fn_perallele_hom)
    outfile.write("HOM Precision (TP / (TP+FP)): {}\n".format(precision_HOM))
    outfile.write("HOM Recall (TP / (TP+FN)): {}\n".format(recall_HOM))
    outfile.write("HOM F1 ( (2*TP)/(2*TP+FP+FN)): {}\n".format(f1_HOM))

    precision_HET = cnt_tp_perallele_het / (cnt_tp_perallele_het + cnt_fp_perallele_het)
    recall_HET = cnt_tp_perallele_het / (cnt_tp_perallele_het + cnt_fn_perallele_het)
    f1_HET = (2*cnt_tp_perallele_het) / (2 * cnt_tp_perallele_het + cnt_fp_perallele_het + cnt_fn_perallele_het)
    outfile.write("Multi-allelic Precision (TP / (TP+FP)): {}\n".format(precision_HET))
    outfile.write("Multi-allelic Recall (TP / (TP+FN)): {}\n".format(recall_HET))
    outfile.write("Multi-allelic F1 ( (2*TP)/(2*TP+FP+FN)): {}\n".format(f1_HET))

    outfile.write("F1 Macro ( 1/2 (HOM F1+Multi-allelic F1)): {}\n".format((f1_HET + f1_HOM) / 2))

    # Additional TLS
    accuracy = cnt_con_perposition / cnt_comparablePositions
    outfile.write("Total Accuracy (Concordance/Comparable Positions): {}\n".format(accuracy))

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
    accuracy = (cnt_con_perposition + cnt_implicit_agreement) / (
        cnt_comparablePositions + cnt_implicit_agreement
    )
    outfile.write(
        "Global Accuracy ((Concordance+IA)/(Comparable Positions+IA)): {}\n".format(accuracy)
    )

    outfile.write(
        "Global Error Rate (1- Global Accuracy): {}\n".format(1-accuracy)
    )

    #Combined Per Allele Stats
    outfile.write("Per Allele Precision (TP/(TP+FP)): {}\n".format(cnt_tp_perallele/(cnt_tp_perallele+cnt_fp_perallele)))
    outfile.write("Per Allele Recall (TP/(TP+FN)): {}\n".format(cnt_tp_perallele/(cnt_tp_perallele+cnt_fn_perallele)))
    outfile.write("Per Allele F1 (2TP/(2TP+FN+FP)): {}\n".format(2*cnt_tp_perallele/(2*cnt_tp_perallele+cnt_fn_perallele+cnt_fp_perallele)))

    #Multiallelic Classification
    outfile.write("Multi-Allelic TP: {}\n".format(cnt_multiallelic_tp))
    outfile.write("Multi-Allelic TN: {}\n".format(cnt_multiallelic_tn))
    outfile.write("Multi-Allelic FP: {}\n".format(cnt_multiallelic_fp))
    outfile.write("Multi-Allelic FN: {}\n".format(cnt_multiallelic_fn))

    ma_prec = (cnt_multiallelic_tp / (cnt_multiallelic_tp + cnt_multiallelic_fp)) if (cnt_multiallelic_tp + cnt_multiallelic_fp) != 0 else '?'
    outfile.write("Multi-Allelic Precision (TP/(TP+FP)): {}\n".format(ma_prec))
    ma_recall = cnt_multiallelic_tp/(cnt_multiallelic_tp+cnt_multiallelic_fn) if (cnt_multiallelic_tp + cnt_multiallelic_fn) != 0 else '?'
    outfile.write("Multi-Allelic Recall (TP/(TP+FN)): {}\n".format(ma_recall))
    outfile.write("Multi-Allelic F1 (2TP/(2TP+FN+FP)): {}\n".format(2*cnt_multiallelic_tp/(2*cnt_multiallelic_tp+cnt_multiallelic_fp+cnt_multiallelic_fn)))

    outfile.write("Multi-Allelic Specificity (TN/(FP+TN)): {}\n".format(cnt_multiallelic_tn/(cnt_multiallelic_tn+cnt_multiallelic_fp)))
    outfile.write("Multi-Allelic Negative Predictive Value (TN/(FN+TN)): {}\n".format(cnt_multiallelic_tn/(cnt_multiallelic_tn+cnt_multiallelic_fn)))

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
        "referenceContext",
    ],
)

df.to_csv(snakemake.output["dataframe"])
