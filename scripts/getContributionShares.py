import os

import vcfpy
import pandas
import upsetplot

methods = ['medaka','nanopolish','freebayes','gisaid']

#DICT Structure [POSITION] -> [VARIANT] -> [LIST OF METHODS THAT CALLED]
calls = {}

#list is structured such that item #1 is pancov item #2 is freebayes item#3 is medaka and item#4 is nanopolish
iterator = iter(snakemake.input)

#We are interested in:
totalVars = {} #total number of vars called by each method
totalUnique = {} #total number of vars called by each method independent of sample

pancovAddShare = {
    k : 0 for k in methods
} #vars that each method contributes to our method

pancovAddExclusiveShare = {
    k: 0 for k in methods
}


pancovDiscardShare = {
    k: 0 for k in methods
}


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


#get last file

aslist = list(iterator)
ga_path = aslist[-1]
processMethod('gisaid',ga_path)

print(ga_path,len(aslist))

for pc_path in iter(aslist[:-1]):
    fb_path = next(iterator)
    md_path = next(iterator)
    np_path = next(iterator)
    processMethod('pancov',pc_path)
    processMethod('freebayes',fb_path)
    processMethod('medaka',md_path)
    processMethod('nanopolish',np_path)


print('Analyzing files ...')

unionSize = 0

df_data = list()
for position in calls:
    for variant in calls[position]:
        unionSize += 1
        df_data.append(('pancov' in calls[position][variant], 'freebayes' in calls[position][variant],
                        'medaka' in calls[position][variant], 'nanopolish' in calls[position][variant],
                        'gisaid' in calls[position][variant],
                        position, variant))
        
        if 'pancov' in calls[position][variant]: #we use this variant
            mask = (
                1 if 'freebayes' in calls[position][variant] else 0,
                1 if 'medaka' in calls[position][variant] else 0,
                1 if 'nanopolish' in calls[position][variant] else 0,
                1 if 'gisaid' in calls[position][variant] else 0
            )
            if mask == (0,0,1,0):
                pancovAddExclusiveShare['nanopolish'] += 1
            elif mask == (0,1,0,0):
                pancovAddExclusiveShare['medaka'] += 1
            elif mask == (1,0,0,0):
                pancovAddExclusiveShare['freebayes'] += 1
            elif mask == (1, 0, 0,0):
                pancovAddExclusiveShare['gisaid'] += 1
            else:
                for method in methods:
                    if method in calls[position][variant]:
                        pancovAddShare[method] += 1
        else:
            if 'freebayes' in calls[position][variant]:
                pancovDiscardShare['freebayes'] += 1
            if 'medaka' in calls[position][variant]:
                pancovDiscardShare['medaka'] += 1
            if 'nanopolish' in calls[position][variant]:
                pancovDiscardShare['nanopolish'] += 1
            if 'gisaid' in calls[position][variant]:
                pancovDiscardShare['gisaid'] += 1

df = pandas.DataFrame(df_data, columns=('pancov', 'freebayes', 'medaka', 'nanopolish','gisaid', 'position', 'variant'))

only_pancov = df[(df['pancov'] == True) & (df['freebayes'] == False) & (df['medaka'] == False) & (df['nanopolish'] == False) & (df['gisaid'] == False)]
only_pancov.to_csv("only_pancov.csv")

df = df.groupby(['pancov', 'freebayes', 'medaka', 'nanopolish','gisaid']).count()["position"]

upsetplot.plot(df)

import matplotlib
matplotlib.pyplot.savefig(snakemake.output[1])



with open(snakemake.output[0],'w') as outfile:
    for method in totalVars:
        outfile.write('{} detected {} variants (across all samples) \n'.format(method,totalVars[method]))
    for method in totalVars:
        outfile.write('{} detected {} unique variants \n'.format(method,sum(totalUnique[method].values())))
    outfile.write('{} variants were used for the pangenome graph\n'.format(unionSize))
    for method in pancovAddShare:
        outfile.write('{} of the pancov variants originated from {} \n'.format(pancovAddShare[method],method))
    for method in pancovAddExclusiveShare:
        outfile.write('{} of the pancov variants were exclusively found by {} \n'.format(pancovAddExclusiveShare[method],method))
    for method in pancovDiscardShare:
        outfile.write('{} of the {} variants were deemed unworthy by pancov \n'.format(pancovDiscardShare[method],method))
