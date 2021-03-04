import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import parsePileupStrandAwareLight,getCoverage
from Bio import SeqIO
import vcfpy
import pandas as pd

### Skript to perform a vcf based comparison between variant calls


#Read required input
pancovVCF = vcfpy.Reader.from_path(snakemake.input['pancovVCF'])
ivarPseudoVCF = pd.read_csv(snakemake.input['iVarVCF'],sep='\t')
illuminapileup = parsePileupStrandAwareLight(snakemake.input["illuminaPileup"])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input["nanoporePileup"])
reference = SeqIO.read(snakemake.input['ref'],'fasta')

### Step 1: Determine all relevant positions
relevantPositions = set()

variantsPancov = {}
variantsIvar = {}

for record in pancovVCF:
	relevantPositions.add(record.POS)
	variantsPancov[record.POS] = record

for position in ivarPseudoVCF['POS'].unique():
	relevantPositions.add(position)
	variantsIvar[position] = ivarPseudoVCF[ivarPseudoVCF.POS == position]

### Step 2: Process

#Counter Variables for Top Level Stats
cnt_realVariants = 0
cnt_realHETVariants = 0
cnt_concordance = 0
cnt_falsePositives = 0
cnt_falseNegatives = 0
cnt_discordance = 0
cnt_detectedVariants = 0
cnt_comparablePositions = 0
cnt_unscoredPositions = 0
cnt_illuminaDropouts = 0
cnt_nanoporeDropouts = 0

for position in relevantPositions:
	#Determine Nanopore and Illumina Coverage
	nanoporeCoverage = getCoverage(nanoporepileup[position]) if position in nanoporepileup else 0
	illuminaCoverage = getCoverage(illuminapileup[position]) if position in illuminapileup else 0
	#Decide whether the position counts or not
	nanoporeDropout = nanoporeCoverage < config['consensusMinCov']
	illuminaDropout = illuminaCoverage < config['consensusMinCov']


	#Dropouts
	if illuminaDropout or nanoporeDropout:
		cnt_unscoredPositions += 1
		if illuminaDropout:
			cnt_illuminaDropouts += 1
		if nanoporeDropout:
			cnt_nanoporeDropouts += 1
	else:
		pass

	outfile.write('{}\t{}\t{}\t{}\n'.format(
		position,
		reference[int(position-1)], #SeqIO is 0-based
		variantsIvar[position] if position in variantsIvar else 'No Variant calls',
		variantsPancov[position] if position in variantsPancov else 'No Variant calls',
	))
