import json
from shared import *
from Bio import SeqIO

#todo: this is (in parts) insanely hackfixy and unelegant

def parseName(filename):
    #data/auxiliary/pileupAnalysis/nanopolish/corona_38_66/23.pileupanalysis.txt
    split1 = filename.split('/')
    barcode = split1[-1].split('.')[0]
    method = split1[-3]
    run = split1[-2]
    return run,method,barcode

homVariantPositions = json.load(open(snakemake.input['homPositions'],'r'))


#Let's create a dictionary where we track what samples have a variant
variantPositionsRaw = json.load(open(snakemake.input['variantPositions'],'r'))
#we simplify a little bit by ignoring diverse variants
variantPositions = {}
for position in variantPositionsRaw:
    idx = int(position)
    variantPositions[idx] = []
    for variant in variantPositionsRaw[position]:
        for sample in variantPositionsRaw[position][variant]:
            variantPositions[idx].append(sample.split(','))


samplesPerPosition = {}
output = {}

reference = SeqIO.read(snakemake.input['reference'],'fasta')


pileupAnalysisDict = {}
for pileup in snakemake.input['pileups']:
    pileupAnalysisDict[pileup] = parsePileup(pileup)

#process each position where we have a hom variant
for position in homVariantPositions:
    idx = int(position)
    output[idx] = {}
    samplesPerPosition[idx] = 0

    referenceBase = reference[idx-1] #1-based vcf files?

    for pileup in snakemake.input['pileups']:
        #we check if this is a reference bearing sample
        run,method,barcode = parseName(pileup)
        #another unelegant hackfix
        if idx not in variantPositions:
            variantPositions[idx] = []
        for sample in variantPositions[idx]:
            if sample[0] == run and sample[1] == method and sample[2] == barcode:
                #In this case we have a sample that does not bear the reference at the given position
                break
        else:
            pileupAnalysis = pileupAnalysisDict[pileup]
            #print(pileupAnalysis)

            #check for dropouts
            if not idx in pileupAnalysis:
                continue
            
            #check for hom criteria 1&2
            totalReads = sum(pileupAnalysis[idx].values())
            if totalReads < snakemake.config['thresholdHomCall_coverage']:
                continue

            try:
                if pileupAnalysis[idx][referenceBase]/totalReads < snakemake.config['thresholdHomCall']:
                    continue
            except:
                print(idx)
                print(pileupAnalysis[idx])
                raise Exception('bla')




            samplesPerPosition[idx] += 1

            for k,v in pileupAnalysis[idx].items():
                if not k in output[idx]:
                    output[idx][k] = 0
                #print(v)
                output[idx][k] += v

'''
#normalization
normalizedOutput = {}
noReferencePositions = []

for position in homVariantPositions:
    idx = int(position)
    sumOfCounts = sum(output[idx][k] for k in output[idx])
    
    normalizedOutput[idx] = {}
    for k in output[idx]:
        normalizedOutput[idx][k] = output[idx][k]/sumOfCounts
    normalizedOutput[idx]['total'] = sum(normalizedOutput[idx][k] for k in normalizedOutput[idx])  
    normalizedOutput[idx]['samplesPerPosition'] = samplesPerPosition[idx]    
'''

with open(snakemake.output['referenceEmissions'],'w') as outfile:
    json.dump(output,outfile)

