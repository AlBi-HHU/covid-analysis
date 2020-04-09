from pysam import VariantFile
from shared import *
import json

callLabels = {}

vcfFile = VariantFile(snakemake.input['vcfFile'])

pileupAnalysis = {}


#Parse pileup analysis file to a pythonesque structure
with open(snakemake.input['pileupAnalysis'],'r') as infile:
	for line in infile.read().splitlines():
		columns = line.split()
		pos = int(columns[0])
		ref = columns[1]
		spread = columns[2]
		spreadDict = {}
		for val in spread.split(';'):
			entry = val.split('=')
			spreadDict[entry[0]] = int(entry[1])
		pileupAnalysis[pos] = spreadDict

threshold = float(snakemake.config['thresholdHomCall'])

#Decision Making:
#If MEDAKA -> just check if HOM or HET
#If Nanopolish -> Check Coverage at position, use 85% (tune in config) threshold for HOM

for rec in vcfFile.fetch():
	position = rec.pos
	localSpreadDict = pileupAnalysis[position]
	totalReads = sum(localSpreadDict.values())

	callLabels[position] = {}

	for altAllele in rec.alts:

		label = '???'

		altAlleleCount = 0
		# + strand
		altAlleleCount += localSpreadDict[altAllele] if altAllele in localSpreadDict else 0
		# - strand
		altAlleleCount += localSpreadDict[altAllele.lower()] if altAllele.lower() in localSpreadDict else 0

		freq = altAlleleCount/totalReads

		if freq >= threshold:
			label = 'HOM'
		else:
			label = 'NONHOM'
		
		callLabels[position][altAllele] = (label,freq)


with open(snakemake.output[0],'w') as outfile:
	json.dump(callLabels,outfile)
