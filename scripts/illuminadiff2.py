from shared import *
import sys
illuminapileup = parsePileupStrandAwareLight(snakemake.input['illuminaPileup'])

with open(snakemake.output['diffFile'],'w') as outFile,open(snakemake.input['iVarInfo'],'r') as ivarInfoFile,open(snakemake.input['pancovInfo'], 'r') as pancovInfoFile:
	outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Position','Ref','Alt','Recovered','Comment','Pileup'))

	for l in ivarInfoFile.read().splitlines():
		lineData = l.split()
		position = int(lineData[0])
		reference = lineData[1]
		altallele = lineData[2]

		if position in illuminapileup:
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
				outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele,-1,'Did not pass Alex Perl Filter, we ignore it',illuminapileup[position]))
			else:
				recovered = False
				for l2 in pancovInfoFile.read().splitlines():
					lineData2 = l.split()
					position2 = lineData[0]
					reference2 = lineData[1]
					altallele2 = lineData[2]
					if position2 == position and altallele2 == altallele:
						recovered = True
						break
				outFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele,recovered,'',illuminapileup[position]))
		else:
			print('Position {} not covered by illumina pileup (but called in ivar, this is fishy)'.format(position))
			sys.exit(-1)
