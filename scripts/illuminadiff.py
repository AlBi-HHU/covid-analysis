from shared import *

ambiguityChars = {
	'R' : frozenset(('A', 'G')),
	'Y' : frozenset(('C', 'T')),
	'S' : frozenset(('G', 'C')),
	'W' : frozenset(('A', 'T')),
	'K' : frozenset(('T', 'G')),
	'M' : frozenset(('A', 'C'))
}

illuminapileup = parsePileupStrandAwareLight(snakemake.input['illuminaPileup'])
nanoporepileup = parsePileupStrandAwareLight(snakemake.input['nanoporePileup'])

with open(snakemake.output['diffFile'],'w') as outFile, open(snakemake.input['pancovInfo'],'r') as pancovInfoFile:

	outFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Position','Ref','Alt','Rejected','Comment','Illumina','Nanopore'))

	for l in pancovInfoFile.read().splitlines():
		lineData = l.split()
		position = int(lineData[0])
		reference = lineData[1]
		altallele = lineData[2]
		altallele_unmodified = altallele
		#Inspect the non-reference part of two base ambiguities
		if altallele == 'N':
			continue
		if altallele in ambiguityChars:
			(altallele,) = ambiguityChars[altallele]-{reference}
		comment = ''
		reject = False
		if position in illuminapileup and sum(illuminapileup[position].values()) > snakemake.config['illuminaCoverageCutoff']:
			sb = getStrandBias(illuminapileup[position],altallele)
			comment += 'strand bias for allele {}: {}'.format(altallele,sb)
			cov = getCoverage(illuminapileup[position],altallele)
			comment += ' coverage for allele {}: {}'.format(altallele,cov)
			reject =  (cov < snakemake.config['consensusMinCov']) or (min(1-sb,sb) <snakemake.config['consensusStrandBiais'])
			if reject:
				comment = 'REJECTED '+ comment
		else:
			reject = -1
			comment += 'position not covered by illumina reads (dropout?)'

		outFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(position,reference,altallele_unmodified,reject,comment,
		                                                    illuminapileup[position] if position in illuminapileup else 'no illu pileup for pos',
															nanoporepileup[position] if position in nanoporepileup else 'no nano pileup for pos'

		))