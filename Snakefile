#9 -> clean 8 -> mixed

ks = [11,13,15]
methods = ['medaka','nanopolish']
barcodes = [8,9]

rule all:
	input:
		expand('data/auxiliary/graphs/{method}/{barcode}/{k}.gfa',method=methods,barcode=barcodes,k=ks),
		'data/output/tobigram.svg'

include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
