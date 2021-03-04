import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import parsePileupStrandAwareLight,getTotalCoverage
from Bio import SeqIO
import vcfpy
import pandas as pd

### Skript to perform a vcf based comparison between variant calls


#Read required input
pancovVCF = vcfpy.Reader.from_path(snakemake.input['pancovVCF'])
ivarPseudoVCF = pd.read_csv(snakemake.input['iVarVCF'],sep='\t')
ivarPseudoVCF = ivarPseudoVCF[ivarPseudoVCF.PASS != False] #Filter passed vars only
illuminapileup = parsePileupStrandAwareLight(snakemake.input["illuminaPileup"])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input["nanoporePileup"])
reference = SeqIO.read(snakemake.input['ref'],'fasta')

### Step 1: Determine all relevant positions
relevantPositions = set()

recordsNanopore = {}
recordsIllumina = {}

for record in pancovVCF:
	relevantPositions.add(record.POS)
	recordsNanopore[record.POS] = record

for position in ivarPseudoVCF['POS'].unique():
	relevantPositions.add(position)
	recordsIllumina[position] = ivarPseudoVCF[ivarPseudoVCF.POS == position]

### Step 2: Process

#Counter Variables for Top Level Stats
cnt_realVariants = 0#

cnt_realHETVariants = 0
cnt_detectedHETVariants = 0

cnt_realSNP = 0
cnt_detectedSNP = 0

cnt_realINS = 0
cnt_detectedINS = 0

cnt_realDEL = 0
cnt_detectedDEL = 0

cnt_concordance = 0
cnt_falsePositives = 0
cnt_falseNegatives = 0
cnt_discordance = 0
cnt_detectedVariants = 0#
cnt_comparablePositions = 0 #
cnt_unscoredPositions = 0 #
cnt_illuminaDropouts = 0 #
cnt_nanoporeDropouts = 0 #

with open(snakemake.output[0],'w') as outfile:

	fields = ['pos','ref','ivar','nanopore-method']
	outfile.write('\t'.join(fields)+'\n')

	for position in relevantPositions:
		#Determine Nanopore and Illumina Coverage
		nanoporeCoverage = getTotalCoverage(nanoporepileup[position]) if position in nanoporepileup else 0
		illuminaCoverage = getTotalCoverage(illuminapileup[position]) if position in illuminapileup else 0
		#Decide whether the position counts or not
		nanoporeDropout = nanoporeCoverage < snakemake.config['consensusMinCov']
		illuminaDropout = illuminaCoverage < snakemake.config['consensusMinCov']


		#Dropouts
		if illuminaDropout or nanoporeDropout:
			cnt_unscoredPositions += 1
			if illuminaDropout:
				cnt_illuminaDropouts += 1
			if nanoporeDropout:
				cnt_nanoporeDropouts += 1
		else:
			pass

		#Resolve Type and Value of Variant

		nanoporeType = 'None'  # INS,DEL,SNP
		nanoporeValue = ''  # Length of Del / Alt Allele / Insertion Seq
		if position in recordsNanopore:
			nanoporeType = recordsNanopore[position].ALT[0].type
			if nanoporeType == 'INS':
				nanoporeValue = recordsNanopore[position].ALT[0].value[2:] #Ignore first char as this is REF
			elif nanoporeType == 'DEL':
				nanoporeValue = str(len(recordsNanopore[position].REF)-1) #Ignore first char as this is retained
			elif nanoporeType == 'SNV':
				nanoporeValue = recordsNanopore[position].ALT[0].value

		illuminaType = '?' #INS,DEL,SNP
		illuminaValue = '?' #Length of Del / Alt Allele / Insertion Seq

		outfile.write('{}\t{}\t{}\t{}\n'.format(
			position,
			reference[int(position-1)], #SeqIO is 0-based
			recordsIllumina[position]['ALT'].values[0] if position in recordsIllumina else 'No Variant calls',
			nanoporeType+' '+nanoporeValue
		))
	#Output some stats
	cnt_comparablePositions = len(relevantPositions)-cnt_unscoredPositions
	cnt_detectedVariants = len(recordsNanopore)
	cnt_realVariants = len(recordsIllumina)