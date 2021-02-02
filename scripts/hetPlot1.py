import sys
import altair as alt
import pandas as pd
import vcfpy
from shared import *
import os

tuples = []

cutoffHeterozigosity = 0.1

iterator = iter(snakemake.input)

for pancovF, ivarF,nanoporeF,pileupF in zip(snakemake.input["pancov"], snakemake.input["ivar"],snakemake.input['nanopore'],snakemake.input['illuminaPileups']):


    #We will store all positions that are in any of the vcfs/tables here and later add the illumina pileups
    pileupPositions = {}

    #Pankoff
    reader = vcfpy.Reader.from_path(pancovF)
    for record in reader:

        alleleFrequency = record.INFO["VCOV"] / (
            record.INFO["VCOV"] + record.INFO["RCOV"]
        )
        if cutoffHeterozigosity < alleleFrequency < 1 - cutoffHeterozigosity:
            tuples.append((pancovF, record.POS, "pancov", alleleFrequency))
            pileupPositions[record.POS] = record.ALT

    #Nanopolish
    reader = vcfpy.Reader.from_path(nanoporeF)
    for record in reader:

        alleleFrequency = record.INFO["SupportFraction"]
        if cutoffHeterozigosity < alleleFrequency < 1-cutoffHeterozigosity:
            tuples.append((pancovF, record.POS, "nanopolish", alleleFrequency))
            pileupPositions[record.POS] = record.ALT

    #Ajvar
    ivartable = open(ivarF, "r").read().splitlines()[1:]
    illuminafreq = -1
    for il in ivartable:
        d = il.split()
        illuminafreq = float(d[10])
        if cutoffHeterozigosity < alleleFrequency < 1 - cutoffHeterozigosity:
            tuples.append((pancovF, record.POS, "ivar", illuminafreq))
            pileupPositions[record.POS] = record.ALT

    #Add Pileups
    pileup = parsePileupStrandAwareLight(pileupF)
    for pos in pileupPositions:
        pos = int(pos)
        if pos in pileup:
            af = getAlleleFrequency(pileup[pos],pileupPositions[pos])
            tuples.append((pancovF, record.POS, "illumina", af))

    print('added {} potential het sites for file: {}'.format(len(pileupPositions),pancovF))

df = pd.DataFrame(tuples, columns=["file", "pos", "method", "rvt"])

charts = []
os.makedirs(snakemake.output[0],exist_ok=True)
for f in df["file"].unique():

    outfile = f.split('.')[-2].split('/')
    outfile = outfile[-3]+'_'+outfile[-2]
    charts.append(

        alt.Chart(df[df.file == f],title=f).mark_bar().encode(
            x="method:N",
            y=alt.X("rvt:Q", scale=alt.Scale(domain=[0, 1])),
            color="method:N",
            column="pos:O",
            tooltip=["rvt"],
        ).interactive().save(os.path.join(snakemake.output[0],outfile+'.html'))

    )
