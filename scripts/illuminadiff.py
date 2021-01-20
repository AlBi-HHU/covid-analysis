from shared import *

illuminapileup = parsePileupStrandAwareLight(snakemake.input['illuminaPileup'])

with open(snakemake.output['diffFile'],'w') as outFile, open(snakemake.input['pancovInfo'],'r') as pancovInfoFile:

	outFile.write('{}\t{}\t{}\t{}\t{}\n'.format('Position','Ref','Alt','Verified','Illumina'))

	for l in pancovInfoFile.read().splitlines():
		lineData = l.split()
		position = int(lineData[0])
		reference = lineData[1]
		altallele = lineData[2]

		comment = ''
		reject = False
		if position in illuminapileup:
			sb = getStrandBias(illuminapileup[position],altallele)
			comment += 'strand bias for allele {}: {}'.format(altallele,sb)
			cov = getCoverage(illuminapileup[position],altallele)
			comment += ' coverage for allele {}: {}'.format(altallele,cov)
			reject =  (cov < snakemake.config['consensusMinCov']) or (min(1-sb,sb) <snakemake.config['consensusStrandBiais'])
			if reject:
				comment = 'REJECTED '+ comment
		else:
			comment += 'position not covered by illumina reads (dropout?)'

		outFile.write('{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele,reject,comment))