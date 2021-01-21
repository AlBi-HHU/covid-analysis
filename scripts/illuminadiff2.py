from shared import *
import sys

ambiguityChars = {
	'R' : frozenset(('A', 'G')),
	'Y' : frozenset(('C', 'T')),
	'S' : frozenset(('G', 'C')),
	'W' : frozenset(('A', 'T')),
	'K' : frozenset(('T', 'G')),
	'M' : frozenset(('A', 'C'))
}

illuminapileup = parsePileupStrandAwareLight(snakemake.input['illuminaPileup'])
pancovInfoFile= open(snakemake.input['pancovInfo'], 'r').read().splitlines()

with open(snakemake.output['diffFile'],'w') as outFile,open(snakemake.input['iVarInfo'],'r') as ivarInfoFile:
	outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Position','Ref','Alt','Recovered','Comment','Pileup'))

	for l in ivarInfoFile.read().splitlines():
		lineData = l.split()
		position = int(lineData[0])
		reference = lineData[1]
		altalleles = [lineData[2]]
		altallele_unmodified = altallele

		if altallele_unmodified == 'N':
			continue
		if altallele_unmodified in ambiguityChars:
			altalleles = list(ambiguityChars[altallele]-{reference})

		if position in illuminapileup:
			for altallele in altalleles:
				#Alex perl
				sb = getStrandBias(illuminapileup[position],altallele)
				cov = getCoverage(illuminapileup[position],altallele)
				abs = getMinorStrandAbs(illuminapileup[position],altallele)
				fq = getMinorStrandFrequency(illuminapileup[position],altallele)

				reject = False
				if cov <= 10:
					reject = True
				elif cov <= 20:
					if abs < 5:
						reject = True
				elif cov <= 50:
					if abs < 10 and fq < 0.25:
						reject = True
				elif cov <= 100:
					if abs < 15 and fq < 0.15:
						reject = True
				else:
					if fq < 0.1:
						reject = True
				if reject:
					outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele_unmodified,-1,'Allele Component: {} did not pass Alex Perl Filter, we ignore it'.format(altallele),illuminapileup[position]))
					break
			else:
				recovered = False
				for l2 in pancovInfoFile:
					#print(l,l2)
					lineData2 = l2.split()
					position2 = int(lineData2[0])
					reference2 = lineData2[1]
					altallele2 = lineData2[2]
					if position2 == position and altallele2 == altallele_unmodified:
						recovered = True
						break
				outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele_unmodified,recovered,'',illuminapileup[position]))
		else:
			print('Position {} not covered by illumina pileup (but called in ivar, this is fishy)'.format(position))
			sys.exit(-1)
