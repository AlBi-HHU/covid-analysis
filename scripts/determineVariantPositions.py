from pysam import VariantFile
from shared import *
import json

def getIdentifierFromFileName(filename):
	return ','.join([filename.split('/')[-3],filename.split('/')[-1].split('.')[1],filename.split('/')[-1].split('.')[0][7:]])

positions = {}

for vcfFilePath in snakemake.input['vcfFiles']:
	vcfFile = VariantFile(vcfFilePath)
	for rec in vcfFile.fetch():
		position = rec.pos
		alleles = set(rec.alts)

		if not position in positions:
			positions[position] = {}
		for allele in alleles:

			if allele in positions[position]:
				positions[position][allele].append(getIdentifierFromFileName(vcfFilePath))
			else:
				positions[position][allele] = [getIdentifierFromFileName(vcfFilePath)]

with open(snakemake.output[0],'w') as outfile:
	json.dump(positions,outfile)
