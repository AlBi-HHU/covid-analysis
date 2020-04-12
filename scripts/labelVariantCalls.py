from shared import *
import json

callLabels = {}

vcfFile = []
with open(snakemake.input['vcfFile'],'r') as vcffileRaw:
	for line in vcffileRaw.read().splitlines():
		if line.startswith('#'):
			continue
		vcfFile.append(line.split('\t'))

isMedaka = 'medaka' in snakemake.input['vcfFile']

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
coverageThreshold = int(snakemake.config['thresholdHomCall_coverage'])

for rec in vcfFile:
	position = int(rec[1])
	localSpreadDict = pileupAnalysis[position]
	totalReads = sum(localSpreadDict.values())

	#check if Medaka predicts HET
	medakaHETCall = False
	if isMedaka:
		sampleCol = rec[-1]
		calls = sampleCol.split(':')[0].split('/')
		#print(calls)
		if calls[0] != calls[1]:
			medakaHETCall = True

	callLabels[position] = {}
	alts = rec[4].split(',')
	for altAllele in alts:

		label = '???'

		altAlleleCount = 0
		# + strand
		altAlleleCount += localSpreadDict[altAllele] if altAllele in localSpreadDict else 0
		# - strand
		altAlleleCount += localSpreadDict[altAllele.lower()] if altAllele.lower() in localSpreadDict else 0

		freq = altAlleleCount/totalReads

		#If the frequency is above the threshold, and we have a minimum coverage and medaka did not call it as HET, we label it as HOM
		if freq >= threshold and totalReads >= coverageThreshold and medakaHETCall == False:
			label = 'HOM'
		#If medaka labels as HET we can look into it further (mixture model) TODO: Define criterion for Non-Medaka VCFs/pileups when we label as HET
		elif totalReads >= coverageThreshold and medakaHETCall == True:
			label = 'HET'
		else:
			label = 'NONHOM'
		
		callLabels[position][altAllele] = (label,freq,localSpreadDict)


with open(snakemake.output[0],'w') as outfile:
	json.dump(callLabels,outfile)
