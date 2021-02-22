from shared import *

ambiguityChars = {
    "R": frozenset(("A", "G")),
    "Y": frozenset(("C", "T")),
    "S": frozenset(("G", "C")),
    "W": frozenset(("A", "T")),
    "K": frozenset(("T", "G")),
    "M": frozenset(("A", "C")),
}

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
        altallele_unmodified = altallele

        # Inspect the non-reference part of two base ambiguities
        if altallele == "N":
            continue
        if altallele in ambiguityChars:
            (altallele,) = ambiguityChars[altallele] - {reference}

        comment = ""

        if position not in nanoporepileup or (
            position in nanoporepileup
            and sum(nanoporepileup[position].values())
            < snakemake.config["nanoporeCoverageCutoff"]
        ):
            reject = "Nanopore"
            comment += "position not covered by nanopore reads (dropout?)"
        elif position not in illuminapileup or (
            position in illuminapileup
            and sum(illuminapileup[position].values())
            < snakemake.config["illuminaCoverageCutoff"]
        ):
            reject = "Illumina"
            comment += "position not covered by illumina reads (dropout?)"
        else:
            sb = getStrandBias(illuminapileup[position], altallele)
            comment += "strand bias for allele {}: {}".format(altallele, sb)

            cov = getCoverage(illuminapileup[position], altallele)
            comment += " illumina coverage for allele {}: {}".format(altallele, cov)

            if (cov < snakemake.config["consensusMinCov"]) or (
                    min(1 - sb, sb) < snakemake.config["consensusStrandBiais"]
            ):
                reject = "True"
            else:
                reject = "False"

        outFile.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                position,
                reference,
                altallele_unmodified,
                reject,
                comment,
                illuminapileup[position]
                if position in illuminapileup
                else "no illu pileup for pos",
                nanoporepileup[position]
                if position in nanoporepileup
                else "no nano pileup for pos",
            )
        )
