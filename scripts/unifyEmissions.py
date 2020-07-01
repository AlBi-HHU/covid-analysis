import json
from functools import reduce

variantEmissions = json.load(open(snakemake.input['varEmissions'],'r'))
referenceEmissions = json.load(open(snakemake.input['refEmissions'],'r'))

#we can only look at positions where we have both var and ref emissions
positions = set(variantEmissions.keys()).union(set(referenceEmissions.keys()))

unifiedEmissions = {}

for position in positions:
    #determine unified kmer set
    unifiedAlleleSet = reduce( lambda x, y: x.union(y.keys()), [referenceEmissions[position]]+[variantEmissions[position][allele] for allele in variantEmissions[position]], set())
    addition = False

    for allele in unifiedAlleleSet:
        
        if not allele in referenceEmissions[position]:
            addition = True    
            break
    unifiedEmissions[position] = {}
    unifiedEmissions[position]['reference'] = {}


    for allele in unifiedAlleleSet:
        if allele not in referenceEmissions[position]:
            referenceEmissions[position][allele] = 0
        unifiedEmissions[position]['reference'][allele] = referenceEmissions[position][allele]+addition

    for variant in variantEmissions[position]:
        

        addition = False

        for allele in unifiedAlleleSet:
            
            if not allele in variantEmissions[position][variant]:
                addition = True    
                break

        unifiedEmissions[position][variant] = {}

        for allele in unifiedAlleleSet:
            if allele not in variantEmissions[position][variant]:
                variantEmissions[position][variant][allele] = 0
            unifiedEmissions[position][variant][allele] = variantEmissions[position][variant][allele]+addition


    #normalize
    for k in unifiedEmissions[position]:
        total = sum(unifiedEmissions[position][k].values())
        for allele in unifiedEmissions[position][k]:
            unifiedEmissions[position][k][allele] /= total

with open(snakemake.output[0],'w') as outfile:
    json.dump(unifiedEmissions,outfile)
