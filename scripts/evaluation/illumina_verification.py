import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import *


illuminapileup = parsePileupStrandAwareLight(snakemake.input["illuminaPileup"])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input["nanoporePileup"])

with open(snakemake.output["diffFile"], "w") as outFile, open(
    snakemake.input["pancovInfo"], "r"
) as pancovInfoFile:

    outFile.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "Position",
            "Ref",
            "Alt",
            "Rejected",
            "Comment",
            "Illumina",
            "Nanopore",
        )
    )

    for l in pancovInfoFile.read().splitlines():
        lineData = l.split()
        position = int(lineData[0])
        reference = lineData[1]
        altallele = lineData[2]

        # Inspect the non-reference part of two base ambiguities
        if altallele == "N":
            continue

        comment = ""
        status = ""

        if position not in illuminapileup or (
            position in illuminapileup
            and sum(illuminapileup[position].values())< snakemake.config["illuminaCoverageCutoff"]
        ):
            status = "IlluminaDropout"
            comment += "position not sufficiently covered by illumina reads (dropout?)"
        else:
            #Check EACH component for ambiguity:
            checklist = ambiguityLetters_inverted[altallele] if altallele in ambiguityLetters_inverted else [altallele]
            for allele in checklist:
                sb = getStrandBias(illuminapileup[position], allele)
                cov = getCoverage(illuminapileup[position], allele)

                if (cov < snakemake.config["consensusMinCov"]):
                    status = "Rejected"
                    comment = " illumina coverage for component {} too low: {}".format(allele, cov)
                    break
                elif (
                        min(1 - sb, sb) < snakemake.config["pangenomeStrandBiais"]
                ):
                    status = "Rejected"
                    comment = "strand bias for component {} too extreme: {}".format(allele, sb)
                    break
                else:
                    status = "Verified"

        outFile.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                position,
                reference,
                altallele,
                status,
                comment,
                illuminapileup[position]
                if position in illuminapileup
                else "no illu pileup for pos",
                nanoporepileup[position]
                if position in nanoporepileup
                else "no nano pileup for pos",
            )
        )
