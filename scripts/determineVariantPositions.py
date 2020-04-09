from pysam import VariantFile
from shared import *
import json

positions = {}

for vcfFilePath in snakemake.input['vcfFiles']:
	vcfFile = VariantFile(vcfFilePath)
	for rec in vcfFile.fetch():
		position = rec.pos
		alleles = set(rec.alleles)

		if position in positions:
			positions[position] = list(set(positions[position]).union(alleles))
		else:
			positions[position] = list(alleles)

with open(snakemake.output[0],'w') as outfile:
	json.dump(positions,outfile)
