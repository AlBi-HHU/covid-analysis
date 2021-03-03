import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure

from shared import *


illuminapileup = parsePileupStrandAwareLight(snakemake.input["illuminaPileup"])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input["nanoporePileup"])

pancovInfoFile = open(snakemake.input["pancovInfo"], "r").read().splitlines()



def determineRecoveryStatus(position,altallele):
    if position in illuminapileup:

        # Alex perl #TODO: Move elsewhere
        sb = getStrandBias(illuminapileup[position], altallele)
        cov = getCoverage(illuminapileup[position], altallele)
        abs = getMinorStrandAbs(illuminapileup[position], altallele)
        fq = getMinorStrandFrequency(illuminapileup[position], altallele)

        variantFiltered = False
        if cov <= 10:
            variantFiltered = True
        elif cov <= 20:
            if abs < 5:
                variantFiltered = True
        elif cov <= 50:
            if abs < 10 and fq < 0.25:
                variantFiltered = True
        elif cov <= 100:
            if abs < 15 and fq < 0.15:
                variantFiltered = True
        else:
            if fq < 0.1:
                variantFiltered = True

        if variantFiltered:
            return "Filtered","Filtered by Alex SB Filter"

        if position in nanoporepileup:
            if sum(nanoporepileup[position].values()) < snakemake.config["nanoporeCoverageCutoff"]:
                return "NanoporeDropout","Below Threshold {}<{}".format(nanoporepileup[position].values(),snakemake.config["nanoporeCoverageCutoff"])
        else:
            return "NanoporeDropout","Full Dropout"

        for l2 in pancovInfoFile:
            lineData2 = l2.split()
            position2 = int(lineData2[0])
            altallele2 = lineData2[2]
            if position2 == position:
                if altallele2 == altallele:
                    return "Recovered",""
                else:
                    return "Disagreement","Method called: {}".format(altallele2)

        return "Missed","Missed completely"




with open(snakemake.output["diffFile"], "w") as outFile, open(
    snakemake.input["iVarInfo"], "r"
) as ivarInfoFile:
    outFile.write(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "Position",
            "Ref",
            "Alt",
            "Recovered",
            "Comment",
            "Pileup Illu",
            "Pileup Nano",
        )
    )

    for l in ivarInfoFile.read().splitlines():
        lineData = l.split()
        position = int(lineData[0])
        reference = lineData[1]
        altallele = lineData[2]

        if altallele == "N":
            continue

        recovered,comment = determineRecoveryStatus(position,altallele)


        outFile.write(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                position,
                reference,
                altallele,
                recovered, #the status of the allele
                comment, #additional info
                illuminapileup[position],
                nanoporepileup[position] if position in nanoporepileup else "No nanopore pileup",
            )
        )