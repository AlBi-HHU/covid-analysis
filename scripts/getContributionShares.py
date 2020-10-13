import vcfpy
import os

#DICT Structure [POSITION] -> [VARIANT] -> [LIST OF METHODS THAT CALLED]
calls = {}

#list is structured such that item #1 is pancov item #2 is freebayes item#3 is medaka and item#4 is nanopolish
iterator = iter(snakemake.input)

#We are interested in:
totalVars = {} #total number of vars called by each method
totalUnique = {} #total number of vars called by each method independent of sample

pancovAddShare = {} #vars that each method contributes to our method
pancovAddShare['medaka'] = 0
pancovAddShare['nanopolish'] = 0
pancovAddShare['freebayes'] = 0
pancovAddExclusiveShare = {}

#init here for convenience
pancovAddExclusiveShare['medaka'] = 0
pancovAddExclusiveShare['nanopolish'] = 0
pancovAddExclusiveShare['freebayes'] = 0

pancovDiscardShare = {}
pancovDiscardShare['medaka'] = 0
pancovDiscardShare['nanopolish'] = 0
pancovDiscardShare['freebayes'] = 0

def processMethod(method,file):
    vcf_reader = vcfpy.Reader(open(file,'r'))
    for record in vcf_reader:
        if not record.POS in calls:
            calls[record.POS] = {}
        for allele in record.ALT:
            if not allele in calls[record.POS]:
                calls[record.POS][allele] = set()
            calls[record.POS][allele].add(method)
            #Add totalVars counter
            if not method in totalVars:
                totalVars[method] = 0
            totalVars[method] += 1
            #Unique vars per position
            if not method in totalUnique:
                totalUnique[method] = {}
            if not record.POS in totalUnique:
                totalUnique[method][record.POS] = 1

print('Reading files ...')

for pc_path in iterator:
    fb_path = next(iterator)
    md_path = next(iterator)
    np_path = next(iterator)
    processMethod('pancov',pc_path)
    processMethod('freebayes',fb_path)
    processMethod('medaka',md_path)
    processMethod('nanopolish',np_path)

print('Analyzing files ...')

for position in calls:
    for variant in calls[position]:
        if 'pancov' in calls[position][variant]: #we use this variant
            mask = (
                1 if 'freebayes' in calls[position][variant] else 0,
                1 if 'medaka' in calls[position][variant] else 0,
                1 if 'nanopolish' in calls[position][variant] else 0
            )
            if mask == (0,0,1):
                pancovAddExclusiveShare['nanopolish'] += 1
            elif mask == (0,1,0):
                pancovAddExclusiveShare['medaka'] += 1
            elif mask == (1,0,0):
                pancovAddExclusiveShare['freebayes'] += 1
            elif mask == (1, 1, 0):
                pancovAddShare['medaka'] += 1
                pancovAddShare['freebayes'] += 1
            elif mask == (0, 1, 1):
                pancovAddShare['medaka'] += 1
                pancovAddShare['nanopolish'] += 1
            elif mask == (1, 0, 1):
                pancovAddShare['nanopolish'] += 1
                pancovAddShare['freebayes'] += 1
            elif mask == (1,1,1):
                pancovAddShare['nanopolish'] += 1
                pancovAddShare['freebayes'] += 1
                pancovAddShare['medaka'] += 1
        else:
            if 'freebayes' in calls[position][variant]:
                pancovDiscardShare['freebayes'] += 1
            if 'medaka' in calls[position][variant]:
                pancovDiscardShare['medaka'] += 1
            if 'nanopolish' in calls[position][variant]:
                pancovDiscardShare['nanopolish'] += 1

with open(snakemake.output[0],'w') as outfile:
    for method in totalVars:
        outfile.write('{} detected {} variants (across all samples) \n'.format(method,totalVars[method]))
    for method in totalVars:
        outfile.write('{} detected {} unique variants \n'.format(method,sum(totalUnique[method].values())))
    for method in pancovAddShare:
        outfile.write('{} of the pancov variants originated from {} \n'.format(pancovAddShare[method],method))
    for method in pancovAddExclusiveShare:
        outfile.write('{} of the pancov variants were exclusively found by {} \n'.format(pancovAddExclusiveShare[method],method))
    for method in pancovDiscardShare:
        outfile.write('{} of the {} variants were deemed unworthy by pancov \n'.format(pancovDiscardShare[method],method))
