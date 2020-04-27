#9 -> clean 8 -> mixed

configfile: "config.yaml"

ks = config['ks']

methods = ['medaka','nanopolish']
#Uses the IGV sessions which is completely arbitrary, could use any other input file here to get the barcode ids

runs = config['runs'] if config['useSubsetOfRuns'] else glob_wildcards('data/input/{run}/result_hac').run

barcodes = {}

if config['useSubsetOfBarcodes']:
	barcodes = config['barcodes']
else:
	for run in runs:
		barcodes[run] = glob_wildcards('data/input/'+run+'/result_hac/barcode{barcode}.medaka.primertrimmed.vcf').barcode


def getInput():
	inputList = []


	inputList.append([
		'data/output/tobigram.svg',
		'data/auxiliary/interestingHOMPositions.json',
		'data/auxiliary/interestingPositions.json',
		'data/auxiliary/sampleClassification.tsv'
	])


	for run in runs:
		if config['generateGFAs']:		
			inputList += expand('data/auxiliary/graphs/{method}/'+run+'/{barcode}/{k}.gfa',method=methods,barcode=barcodes[run],k=ks)
		if config['performQualityControl']:
			inputList += expand('data/output/softClippedSeqs/{method}/'+run+'/{barcode}.html',method=methods,barcode=barcodes[run])
		if config['generateKmerProfiles']:		
			inputList += expand('data/output/kmerHistograms/{method}/'+run+'/{barcode}_{k}.svg',method=methods,barcode=barcodes[run],k=ks)
		if config['performCorrections']:
                        inputList += expand('data/auxiliary/{corrections_type}/{method}/'+run+'/{k}/{barcode}.fasta',method=methods,barcode=barcodes[run],k=config["kmerUseForCorrections"][-1], corrections_type=["corrections", "corrections_clip"])
		if config['generatePangenome']:
                        inputList += ["data/auxiliary/pangenome/all_reads_and_ref.paths.gfa"]
		inputList += expand('data/output/IgvSessions/{method}/'+run+'/{barcode}.igv.xml',method=methods,barcode=barcodes[run])

	return inputList

rule all:
	input:
		getInput()



include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
include: 'rules/variantAnalysis.snk'
include: 'rules/pangenome.snk'
