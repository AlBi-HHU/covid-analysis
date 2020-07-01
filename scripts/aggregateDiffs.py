import json

output = {}

for f in snakemake.input:
    for l in open(f,'r').read().splitlines()[1:]:
        split = l.split()
        idx1 = split[0]
        idx2 = split[1]+'->'+split[2]
        if not idx1 in output:
            output[idx1] = {}
        if not idx2 in output[idx1]:
            output[idx1][idx2] = []
        output[idx1][idx2].append(f)

json.dump(output,open(snakemake.output[0],'w'))
