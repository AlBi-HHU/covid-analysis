#9 -> clean 8 -> mixed

configfile: "config.yaml"

ks = config['ks']

methods = ['medaka','nanopolish']
#Uses the IGV sessions which is completely arbitrary, could use any other input file here to get the barcode ids
barcodes= glob_wildcards('data/input/result_hac/barcode{barcode}.igv.xml').barcode


rule all:
	input:
		expand('data/auxiliary/graphs/{method}/{barcode}/{k}.gfa',method=methods,barcode=barcodes,k=ks),
		expand('data/output/softClippedSeqs/{method}/{barcode}.html',method=methods,barcode=barcodes),
		expand('data/output/kmerHistograms/{method}/{barcode}_{k}.svg',method=methods,barcode=barcodes,k=ks),
		'data/output/tobigram.svg',
		'data/auxiliary/interestingHOMPositions.json',
		'data/auxiliary/interestingPositions.json'
include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
include: 'rules/variantAnalysis.snk'
