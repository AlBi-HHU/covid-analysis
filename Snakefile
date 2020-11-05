import sys
import json

configfile: "config.yaml"

ks = config['ks']

methods = ['medaka']
#Uses the IGV sessions which is completely arbitrary, could use any other input file here to get the barcode ids

runs = config['runs'] if config['useSubsetOfRuns'] else glob_wildcards('data/input/{run}').run

#describes the ending of the vcf files as output by the ARTIC pipeline
vcf_suffix = config['vcf_suffix']
bam_suffix = config['bam_suffix']

barcodes = {}

if config['useSubsetOfBarcodes']:
    barcodes = config['barcodes']
    print("Using the following subset of barcodes: {}".format(barcodes))
else:
    for run in runs:
        barcodes[run] = glob_wildcards('data/input/'+run+'/barcode{barcode}.medaka.'+vcf_suffix).barcode

NWids = glob_wildcards('data/input/gisaidseqs/Germany_NW-HHU-{id}.fasta').id


### REMOVE LATER
gisaidMapping = {}
gisaidMappingInverse = {}

with open('data/input/mappingRunsGisaid.csv', 'r') as infile:
    lines = infile.read().splitlines()
    for l in lines:
        data = l.split()
        run = data[0]
        barcode = int(data[1])
        file = ''
        gisaidID = '-1'
        if len(data) != 2: #entry in the table
            file = data[2]
            gisaidID = file.split('_')[1].split('-')[-1]
        if not run in gisaidMapping:
            gisaidMapping[run] = {}
        gisaidMapping[run][barcode] = file
        gisaidMappingInverse[gisaidID] = run+'_'+barcode

def getGisaidFile(run,barcode):
    if run in gisaidMapping:
        if barcode in gisaidMapping[run]:
            return gisaidMapping[run][barcode]
    return None

for run in runs:
    for barcode in barcodes[run]:
        gisaidpath = getGisaidFile(run,int(barcode))
        print(gisaidpath)
        '''
        if gisaidpath:
            if os.path.isfile(os.path.join('data/input/gisaidseqs',gisaidpath)):
                continue
            else:
                print('No sequence for gisaid id: {} (run: {} bc: {})'.format(gisaidpath,run,barcode))
        else:
            print('No mapping for run: {} bc: {}'.format(run,barcode))
        '''
### REMOVE LATER END

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

        #if config['generateKmerProfiles']:
        #    inputList += expand('data/output/kmerHistograms/{method}/'+run+'/{barcode}_{k}.svg',method=methods,barcode=barcodes[run],k=ks)

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
            #Realignments (removed for now because this is misleading)
            #inputList += expand('data/output/IgvSessions/realignment/{method}/'+run+'/{k}/{barcode}.igv.xml', method=methods, k=ks, barcode=barcodes[run])
            
            if config['discovery_fb']:
                inputList += expand('data/auxiliary/discovery/freebayes/{method}/'+run+'/{k}/{barcode}.vcf', method=methods, k=ks, barcode=barcodes[run])
                #inputList += expand('data/output/discovery/annotatedAggregatedDiffs_freebayes_{k}.json',k=ks)
            #if config['discovery_ctx']:
                #inputList += expand('data/auxiliary/discovery/cortex/{method}/'+run+'/{k}/{barcode}.vcf',method=methods,k=ks,barcode=barcodes[run])
                #inputList += expand('data/output/discovery/annotatedAggregatedDiffs_cortex_{k}.json',k=ks)

        if config['consensus']:
            inputList += expand('data/output/consensus/{method}/'+run+'/{barcode}/consensus.fasta', method=methods, k=ks, barcode=barcodes[run])
            #inputList += expand('data/output/consensus/{method}/curation.xlsx' method=methods)

        if config['VarAnnotSnpEff']:
            inputList += expand('data/auxiliary/pangenome_vc/{method}/'+run+'/{barcode}/filter.annoted.vcf', method=methods, barcode=barcodes[run])

        '''
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
        '''

        #Always create the igv sessions for our input
        inputList += expand('data/output/IgvSessions/{method}/'+run+'/{barcode}.igv.xml',method=methods,barcode=barcodes[run])
		

        #Debug/Eval Stuff
        #inputList += ['data/auxiliary/pangenome_vc/contrib.txt']
        inputList += ['data/output/evaluation/comparisonFastaBased/nanopolish.eval']
        inputList += ['data/output/evaluation/comparisonFastaBased/medaka.eval']
        inputList += ['data/output/evaluation/comparisonFastaBased/manualCuration.eval']

        for id in NWids:
            inputList += ['data/auxiliary/evaluation/consensusVariantExtraction/gisaid/'+id+'.info']

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
include: 'rules/pangenome.snk'
include: 'rules/discovery_shared.snk'
include: 'rules/discovery_fb.snk'
include: 'rules/pangenome_variant_call.snk'
include: 'rules/pangenome_eval.snk'

#include: 'rules/discovery_ctx.snk'


