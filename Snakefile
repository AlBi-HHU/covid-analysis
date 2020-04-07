#9 -> clean 8 -> mixed

ks = [11,13,15,17]
methods = ['medaka','nanopolish']
barcodes = ['08','09']

rule all:
	input:
		expand('data/auxiliary/graphs/{method}/{barcode}/{k}.gfa',method=methods,barcode=barcodes,k=ks),
		expand('data/output/softClippedSeqs/{method}/{barcode}.html',method=methods,barcode=barcodes),
		expand('data/output/kmerHistograms/{method}/{barcode}_{k}.svg',method=methods,barcode=barcodes,k=ks),
		'data/output/tobigram.svg'

include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
