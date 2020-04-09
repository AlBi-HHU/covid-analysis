from pysam import VariantFile
from shared import *
import json

callLabels = {}

vcfFile = VariantFile(snakemake.input['vcfFile'])
threshold = snakemake.config['thresholdHomCall']

#Decision Making:
#If MEDAKA -> just check if HOM or HET
#If Nanopolish -> Check Coverage at position, use 85% threshold for HOM


for rec in vcfFile.fetch():
	position = rec.pos


with open(snakemake.output[0],'w') as outfile:
	json.dump(callLabels,outfile)
