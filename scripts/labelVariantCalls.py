from pysam import VariantFile
from shared import *
import json

callLabels = {}

vcfFile = VariantFile(snakemake.input['vcfFile'])


# We take three sets of information into account: 
# 


for rec in vcfFile.fetch():
	position = rec.pos


with open(snakemake.output[0],'w') as outfile:
	json.dump(callLabels,outfile)
