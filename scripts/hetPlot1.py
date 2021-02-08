import sys
import altair as alt
import pandas as pd
import vcfpy
from shared import *
import os
from functools import reduce

tuples = []
difftuples = []
cutoffHeterozigosity = 0.1

iterator = iter(snakemake.input)

for pancovF, ivarF, nanoporeF, pileupF in zip(
    snakemake.input["pancov"],
    snakemake.input["ivar"],
    snakemake.input["nanopore"],
    snakemake.input["illuminaPileups"],
):

    localtuples = []

    # We will store all positions that are in any of the vcfs/tables here and later add the illumina pileups
    pileupPositions = {}

    # Pankoff
    reader = vcfpy.Reader.from_path(pancovF)
    for record in reader:

        alleleFrequency = record.INFO["VCOV"] / (
            record.INFO["VCOV"] + record.INFO["RCOV"]
        )
        localtuples.append(
            (pancovF, record.POS, "pancov", alleleFrequency, record.ALT[0].value)
        )
        pileupPositions[record.POS] = record.ALT[0].value

    # Nanopolish
    reader = vcfpy.Reader.from_path(nanoporeF)
    for record in reader:

        alleleFrequency = record.INFO["SupportFraction"]
        localtuples.append(
            (pancovF, record.POS, "nanopolish", alleleFrequency, record.ALT[0].value)
        )
        pileupPositions[record.POS] = record.ALT[0].value

    # Ajvar
    ivartable = open(ivarF, "r").read().splitlines()[1:]
    illuminafreq = -1
    for il in ivartable:
        d = il.split()
        pos = d[1]
        altallele = d[3]
        illuminafreq = float(d[10])
        localtuples.append((pancovF, pos, "ivar", illuminafreq, altallele))
        pileupPositions[pos] = altallele

    # Add Pileups
    pileup = parsePileupStrandAwareLight(pileupF)
    for pos in pileupPositions:
        ipos = int(pos)
        if ipos in pileup:
            af = getAlleleFrequency(pileup[ipos], pileupPositions[pos])
            if (
                snakemake.config["thresholdHomCall"]
                < af
                < 1 - snakemake.config["thresholdHomCall"]
            ):
                localtuples.append((pancovF, pos, "illumina", af, "none"))
                for tuple in filter(
                    lambda x: x[1] == pos and x[2] != "illumina", localtuples
                ):
                    difftuples.append(
                        (tuple[0], tuple[2], tuple[3], af)
                    )
            else:
                localtuples = list(filter(lambda x: x[1] != pos, localtuples))
        else:
            pass
            # print(ipos)
    # print('added {} potential het sites for file: {}'.format(len(pileupPositions),pancovF))
    tuples += localtuples

df = pd.DataFrame(tuples, columns=["file", "pos", "method", "alleleFreq", "call"])
diff_df = pd.DataFrame(difftuples, columns=["file", "method", "truth", "alleleFreq"])
diff_df["square_diff"] = (diff_df["truth"] - diff_df["alleleFreq"]).pow(2)
diff_df["abs_diff"] = (diff_df["truth"] - diff_df["alleleFreq"]).abs()

diff_df["diff"] = (diff_df["truth"] - diff_df["alleleFreq"]).abs()

diff_df.to_csv(snakemake.output["resume"])
for (method, data) in diff_df.groupby("method"):
    print(method, data.describe())

os.makedirs(snakemake.output[0], exist_ok=True)
for f in df["file"].unique():

    outfile = f.split(".")[-2].split("/")
    outfile = outfile[-3] + "_" + outfile[-2]

    alt.Chart(df[df.file == f], title=f).mark_bar().encode(
        x="method:N",
        y=alt.X("alleleFreq:Q", scale=alt.Scale(domain=[0, 1])),
        color="method:N",
        column="pos:O",
        tooltip=["alleleFreq", "call"],
    ).save(os.path.join(snakemake.output["hetfolder"], outfile + ".html"))

plots = []
for m in ["pancov", "ivar", "nanopolish"]:
    base = alt.Chart(diff_df[diff_df.method == m], title=m)

    bar = base.mark_bar().encode(
        x=alt.X("diff:Q", bin=True, scale=alt.Scale(domain=[0, 1])), y="count()"
    )

    rule = base.mark_rule(color="red").encode(x="mean(diff):Q", size=alt.value(5))

    plots.append(bar + rule)

    base = alt.Chart(diff_df[diff_df.method == m], title=m + "_abs")

    bar = base.mark_bar().encode(
        x=alt.X("diff_abs:Q", bin=True, scale=alt.Scale(domain=[-1, 1])), y="count()"
    )

    rule = base.mark_rule(color="red").encode(x="mean(diff_abs):Q", size=alt.value(5))

    plots.append(bar + rule)

reduce(alt.vconcat, plots).save(snakemake.output["overview"])
# break
