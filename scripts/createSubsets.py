import json
import random

output = {}

totalNumberOfSamples = len(snakemake.input)
for subdivision in list(range(snakemake.config['cohortSubdivisions']))+['all']:
    cohortSize = '?'
    if not subdivision == 'all':
        cohortSize = int(totalNumberOfSamples / snakemake.config['cohortSubdivisions']) * (subdivision+1)
    else:
        cohortSize = totalNumberOfSamples
    output[cohortSize] = []

    if subdivision == 'all':
        output[cohortSize].append(snakemake.input)
    else:
        for _ in range(snakemake.config['subdivisionSampleSets']):
            output[cohortSize].append(random.sample(snakemake.input,cohortSize))

with open(snakemake.output[0],'w') as outfile:
    json.dump(output,outfile)


'''
for cohortSize in output:
    for idx,s in enumerate(output[cohortSize]):
        with open(snakemake.output['path_prefix']+'/'+str(cohortSize)+'/'+str(idx)+'.sampleSet','w') as outfile:
            output.write('\n'.join(s))
'''
