import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import parsePileupStrandAwareLight,getTotalCoverage,ambiguityLetters,ambiguityLetters_inverted,getCoverage,getMinorStrandAbs,getMinorStrandFrequency,alexSBFilter
from Bio import SeqIO
import vcfpy
import pandas as pd

### Skript to perform a vcf based comparison between variant calls


# Counter Variables for Top Level Stats
cnt_realVariants = 0  #

cnt_realHETSNPs = 0  #
cnt_detectedHETSNPs = 0  #

cnt_realSNP = 0  #
cnt_detectedSNP = 0  #

cnt_realINS = 0  #
cnt_detectedINS = 0  #

cnt_realDEL = 0  #
cnt_detectedDEL = 0  #

cnt_concordance = 0  #
cnt_falsePositives = 0
cnt_falseNegatives = 0
cnt_discordance = 0  #
cnt_detectedVariants = 0  #
cnt_relevantPositions = 0 #
cnt_comparablePositions = 0  #
cnt_unscoredPositions = 0  #
cnt_illuminaDropouts = 0  #
cnt_nanoporeDropouts = 0  #


#We keep track of tuples in addition to construct a pandas dataframe for easy transformation into altair plots
dataTuples = []

with open(snakemake.output['text'],'w') as outfile:

	#Read required input
	reference = SeqIO.read(snakemake.input['ref'],'fasta')

	#Input comes in blocks of fours
	data = iter(snakemake.input['comparisonFiles'])
	for ivarPseudoVCF in data:

		fid = ivarPseudoVCF.split('/')[-1]

		#Write down file name to keep output file readable
		outfile.write(fid+'\n')

		#Variable Reassignment to keep code more readable
		ivarPseudoVCF = pd.read_csv(ivarPseudoVCF, sep='\t')
		ivarPseudoVCF = ivarPseudoVCF[ivarPseudoVCF.PASS != False]  # Filter passed vars only

		pancovVCF = vcfpy.Reader.from_path(next(data))

		illuminapileup = parsePileupStrandAwareLight(next(data))
		nanoporepileup = parsePileupStrandAwareLight(next(data))

		### Step 1: Determine all relevant positions
		relevantPositions = set()

		recordsNanopore = {}
		recordsIllumina = {}

		for record in pancovVCF:
			relevantPositions.add(record.POS)
			recordsNanopore[record.POS] = record

		for position in ivarPseudoVCF['POS'].unique():
			relevantPositions.add(position)

			record = ivarPseudoVCF[ivarPseudoVCF.POS == position]
			altallele = record['ALT'].values[0]
			#if we have a deletion or such we ignore it for the sb test
			if altallele.startswith('-') or altallele.startswith('+'):
				recordsIllumina[position] = record
				break

			#otherwise we apply the strand bias filter test
			components = ambiguityLetters_inverted[altallele] if altallele in ambiguityLetters_inverted else [altallele]

			for component in components:

				cov = getCoverage(illuminapileup[position], altallele)
				abs = getMinorStrandAbs(illuminapileup[position], altallele)
				fq = getMinorStrandFrequency(illuminapileup[position], altallele)

				if alexSBFilter(cov, abs, fq):
					#print(altallele,cov,abs,fq)
					break
			else:
				recordsIllumina[position] = record

		### Step 2: Process

		fields = ['pos','ref','ivar','nanopore-method']
		outfile.write('\t'.join(fields)+'\n')

		for position in relevantPositions:

			#Track some properties for easier visualization of results later on


			bool_falseNegative = False
			bool_falsePositive = False
			bool_heterozygousIllu = False
			bool_heterozygousNano = False
			bool_concordance = False


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


			#Resolve Type and Value of Variant

			nanoporeType = 'Not Called'  # INS,DEL,SNP
			nanoporeValue = ''  # Length of Del / Alt Allele / Insertion Seq
			if position in recordsNanopore:
				nanoporeType = recordsNanopore[position].ALT[0].type
				if nanoporeType == 'INS':
					nanoporeValue = recordsNanopore[position].ALT[0].value[2:] #Ignore first char as this is REF
					cnt_detectedINS += 1
				elif nanoporeType == 'DEL':
					nanoporeValue = str(len(recordsNanopore[position].REF)-1) #Ignore first char as this is retained
					cnt_detectedDEL += 1
				elif nanoporeType == 'SNV':
					if nanoporeValue in ambiguityLetters_inverted:
						bool_heterozygousNano = True
						cnt_detectedHETSNPs += 1
					nanoporeValue = recordsNanopore[position].ALT[0].value
					cnt_detectedSNP += 1

			illuminaType = 'Not Called' #INS,DEL,SNP
			illuminaValue = '' #Length of Del / Alt Allele / Insertion Seq

			if position in recordsIllumina:
				altval = recordsIllumina[position]['ALT'].values[0]
				refval = recordsIllumina[position]['REF'].values[0]
				altfreq = float(recordsIllumina[position]['ALT_FREQ'].values[0])
				if altval.startswith('-'):
					illuminaType = 'DEL'
					cnt_realDEL += 1
				elif altval.startswith('+'):
					illuminaType = 'INS'
					cnt_realINS += 1
				else:
					illuminaType = 'SNV'
					cnt_realSNP += 1
				if illuminaType == 'INS':
					illuminaValue = altval[2:] #Don't use the +
				elif illuminaType == 'DEL':
					illuminaValue = str(len(altval)-1) #Ignore first char as this is retained
				elif illuminaType == 'SNV':
					illuminaValue = altval
					if snakemake.config['thresholdHomCall'] <= altfreq <= (1-snakemake.config['thresholdHomCall']):
						bool_heterozygousIllu = True
						cnt_realHETSNPs += 1
						illuminaValue = ambiguityLetters[frozenset((altval,refval))]

			if not(illuminaDropout or nanoporeDropout):
				if (illuminaType == nanoporeType) and (illuminaValue == nanoporeValue):
					bool_concordance = True
					cnt_concordance += 1
				else:
					cnt_discordance += 1
					if (position in recordsIllumina )and (not position in recordsNanopore):
						bool_falseNegative = True
						cnt_falseNegatives += 1
					if (position in recordsNanopore )and (not position in recordsIllumina):
						bool_falsePositive = True
						cnt_falsePositives += 1


			#Write Text Output
			outfile.write('{}\t{}\t{}\t{}\n'.format(
				position,
				reference[int(position-1)], #SeqIO is 0-based
				illuminaType+' '+illuminaValue,
				nanoporeType+' '+nanoporeValue
			))

			#Write to tuples
			dataTuples.append((
				fid,
				position,
				reference[int(position-1)],
				illuminaType,
				illuminaValue,
				nanoporeType,
				nanoporeValue,
				nanoporeDropout,
				illuminaDropout,
				bool_falseNegative,
				bool_falsePositive,
				bool_heterozygousIllu,
				bool_heterozygousNano,
				bool_concordance,
				illuminapileup[position] if position in illuminapileup else 'Dropout',
				nanoporepileup[position] if position in nanoporepileup else 'Dropout'
			))



		#Calculate some additional stats
		cnt_detectedVariants += len(recordsNanopore)
		cnt_realVariants += len(recordsIllumina)
		cnt_relevantPositions += len(relevantPositions)

	cnt_comparablePositions += cnt_relevantPositions-cnt_unscoredPositions

	outfile.write('Real (iVar) Variants:\t{} \n'.format(cnt_realVariants))
	outfile.write('Real (iVar) SNPs:\t{} \n'.format(cnt_realSNP))
	outfile.write('Detected SNPs:\t{} \n'.format(cnt_detectedSNP))
	outfile.write('Real (iVar) HET SNPs:\t{} \n'.format(cnt_realHETSNPs))
	outfile.write('Detected HET SNPs:\t{} \n'.format(cnt_detectedHETSNPs))
	outfile.write('Real (iVar) INS:\t{} \n'.format(cnt_realINS))
	outfile.write('Detected INS:\t{} \n'.format(cnt_detectedINS))
	outfile.write('Real (iVar) DEL:\t{} \n'.format(cnt_realDEL))
	outfile.write('Detected DEL:\t{} \n'.format(cnt_detectedDEL))
	outfile.write('Concordance:\t{}\t of\t{}\t comparable positions \n'.format(cnt_concordance,cnt_comparablePositions))
	outfile.write('FP:\t{} \n'.format(cnt_falsePositives))
	outfile.write('FN:\t{} \n'.format(cnt_falseNegatives))
	outfile.write('Discordance: \t{}\t of \t{}\t comparable positions \n'.format(cnt_discordance,cnt_comparablePositions))
	outfile.write('Detected Variants: \t{} \n'.format(cnt_detectedVariants))
	outfile.write('Unscored Positions: \t{}\t (\t{}\t Illumina and \t{}\t Nanopore Dropouts) \n'.format(cnt_unscoredPositions,cnt_illuminaDropouts,cnt_nanoporeDropouts))


#Write Pandas Dataframe
df = pd.DataFrame(dataTuples,columns=[
	'sample',
	'position',
	'reference',
	'illuminaType',
	'illuminaValue',
	'nanoporeType',
	'nanoporeValue',
	'nanoporeDropout',
	'illuminaDropout',
	'falseNegative',
	'falsePositive',
	'realHET',
	'detectedHET',
	'concordance',
	'illuminaPileup',
	'nanoporePileup'
])

df.to_csv(snakemake.output['dataframe'])