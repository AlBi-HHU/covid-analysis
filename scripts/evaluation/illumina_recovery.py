import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure

from shared import *


illuminapileup = parsePileupStrandAwareLight(snakemake.input["illuminaPileup"])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input["nanoporePileup"])

pancovInfoFile = open(snakemake.input["pancovInfo"], "r").read().splitlines()

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
        altalleles = lineData[2]

        if altallele_unmodified == "N":
            continue

        if position in illuminapileup:
            for altallele in altalleles:
                # Alex perl
                sb = getStrandBias(illuminapileup[position], altallele)
                cov = getCoverage(illuminapileup[position], altallele)
                abs = getMinorStrandAbs(illuminapileup[position], altallele)
                fq = getMinorStrandFrequency(illuminapileup[position], altallele)

                reject = False
                if cov <= 10:
                    reject = True
                elif cov <= 20:
                    if abs < 5:
                        reject = True
                elif cov <= 50:
                    if abs < 10 and fq < 0.25:
                        reject = True
                elif cov <= 100:
                    if abs < 15 and fq < 0.15:
                        reject = True
                else:
                    if fq < 0.1:
                        reject = True

                if reject:
                    outFile.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            position,
                            reference,
                            altalleles,
                            "Reject",
                            "Allele Component: {} did not pass Alex Perl Filter, we ignore it".format(
                                altallele
                            ),
                            illuminapileup[position],
                            nanoporepileup[position]
                            if position in nanoporepileup
                            else "No nanopore pileup",
                        )
                    )
                    break
            else:
                recovered = "False"

                for l2 in pancovInfoFile:
                    # print(l,l2)
                    lineData2 = l2.split()
                    position2 = int(lineData2[0])
                    reference2 = lineData2[1]
                    altallele2 = lineData2[2]
                    if position2 == position and altallele2 == altallele_unmodified:
                        recovered = "True"
                        break

                if position in nanoporepileup and sum(nanoporepileup[position].values()) < snakemake.config["nanoporeCoverageCutoff"]:
                    recovered = "Nanopore"

                outFile.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        position,
                        reference,
                        altalleles,
                        recovered,
                        "",
                        illuminapileup[position],
                        nanoporepileup[position]
                        if position in nanoporepileup
                        else "No nanopore pileup",
                    )
                )
        else:
            print(
                "Position {} not covered by illumina pileup (but called in ivar, this is fishy)".format(
                    position
                )
            )
            sys.exit(-1)
