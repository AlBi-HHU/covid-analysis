import json

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


def getInput(wildcards):
    inputList = [
                      'data/output/tobigram.svg',
                      'data/auxiliary/interestingHOMPositions.json',
                      'data/auxiliary/interestingPositions.json',
                      'data/auxiliary/sampleClassification.tsv'
                      ]

    for run in runs:
        if config['generateGFAs']:
            inputList += expand('data/auxiliary/graphs/{method}/'+run+'/{barcode}/{k}.gfa',method=methods,barcode=barcodes[run],k=ks)

        if config['performQualityControl']:
            inputList += expand('data/output/softClippedSeqs/{method}/'+run+'/{barcode}.html',method=methods,barcode=barcodes[run])

        if config['generateKmerProfiles']:
            inputList += expand('data/output/kmerHistograms/{method}/'+run+'/{barcode}_{k}.svg',method=methods,barcode=barcodes[run],k=ks)

        if config['performCorrections']:
            inputList += expand('data/output/corrections/{basecalling}/'+run+'/{k}/{barcode}.svg',
                                basecalling=methods,
                                barcode=barcodes[run],
                                k=ks
                                )

        if config['generatePangenome']:
           inputList += ["data/auxiliary/pangenome/pangenome.gfa"]

        if config['pangenomeVariantCalling']:
            inputList += expand('data/auxiliary/pangenome_vc/{method}/'+run+'/{barcode}/filter.vcf', method=methods, barcode=barcodes[run])

        if config['discovery']:
            #TODO: Rename to "Discovery"? but is mandatory ...
            #Realignments
            inputList += expand('data/output/IgvSessions/realignment/{method}/'+run+'/{k}/{barcode}.igv.xml', method=methods, k=ks, barcode=barcodes[run])
            
            if config['discovery_fb']:
                inputList += expand('data/auxiliary/discovery/freebayes/{method}/'+run+'/{k}/{barcode}.vcf', method=methods, k=ks, barcode=barcodes[run])
                #inputList += expand('data/output/discovery/annotatedAggregatedDiffs_freebayes_{k}.json',k=ks)
            if config['discovery_ctx']:
                inputList += expand('data/auxiliary/discovery/cortex/{method}/'+run+'/{k}/{barcode}.vcf',method=methods,k=ks,barcode=barcodes[run])
                #inputList += expand('data/output/discovery/annotatedAggregatedDiffs_cortex_{k}.json',k=ks)

        if config['consensus']:
            inputList += expand('data/output/consensus/{method}/'+run+'/{barcode}/consensus.fasta', method=methods, k=ks, barcode=barcodes[run])
            #inputList += expand('data/output/consensus/{method}/curation.xlsx' method=methods)
                
        if config['performMethodEvaluation']:
            sampleSetsFile = checkpoints.createSubsets.get().output
            with open(str(sampleSetsFile),'r') as infile:
                sampleSets = json.load(infile)
                for cohortSize in sampleSets:
                    for idx,sampleSet in enumerate(sampleSets[cohortSize]):
                        #Very ugly in my opinion
                        for f in sampleSet:
                            #working with 'data/auxiliary/corrections/{method}/{run}/21/{barcode}.fasta'
                            split = f.split('/')
                            barcode = split[-1].split('.')[0]
                            run = split[-3]
                            method = split[-4]
                            fstring = 'data/auxiliary/pangenome_vc/'+cohortSize+'/'+str(idx)+'/'+method+'/'+run+'/'+barcode+'/variant.vcf'
                            inputList.append(fstring)
        #Always create the igv sessions for our input
        inputList += expand('data/output/IgvSessions/{method}/'+run+'/{barcode}.igv.xml',method=methods,barcode=barcodes[run])
		

        #Debug/Eval Stuff
        inputList += ['data/auxiliary/pangenome_vc/contrib.txt']

    return inputList


rule all:
    input:
        getInput

include: 'rules/shared.snk'
include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
include: 'rules/variantAnalysis.snk'
include: 'rules/consensus.snk'

if config['performMethodEvaluation']:
    include: 'rules/pangenome_eval.snk'
else:
    include: 'rules/pangenome.snk'
    include: 'rules/discovery_shared.snk'
    include: 'rules/discovery_fb.snk'
    include: 'rules/discovery_ctx.snk'
    include: 'rules/pangenome_variant_call.snk'

