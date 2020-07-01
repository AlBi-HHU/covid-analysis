import json

emissions = json.load(open(snakemake.input['emissions'],'r'))
variantPositions = json.load(open(snakemake.input['variantPositions'],'r'))

output = {}

for p in variantPositions:
    if not p in emissions:
        output[p] = variantPositions[p]

with open(snakemake.output[0],'w') as outfile:
    json.dump(output,outfile)

