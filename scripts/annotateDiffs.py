from shared import *
import json


def parseFilename(f):
    d = f.split('/')
    if d[3] == 'freebayes':
        return (d[-1].split('.')[0],d[-3],d[-4])
    else:
        return (d[-1].split('.')[0],d[-2],d[-3])

diffs = json.load(open(snakemake.input['diffs'],'r'))
for pos in diffs:
    for variant in diffs[pos]:
        print('Annotating {} {} ...'.format(pos,variant))
        #Determine Length of Variant
        tmp = variant.split('->')
        variantLength = -1
        if (tmp[0] == 'Missing'): #We detected  a new variant
            pass
        else:
            continue #We missed an existing variant ... doesn't matter
        variantLength = max(len(tmp[1]),len(tmp[2]))
        strandBiases = []
        coverageValues = []
        for f in snakemake.input['pileups']:
            fp = parseFilename(f)
            for fn in diffs[pos][variant]:
                fnp = parseFilename(fn)
                if fnp == fp:
                    #get strand bias and coverage across the entire length of variant
                    strandBias = 0
                    coverage = 0
                    posAsInt = int(pos)
                    for p in range(posAsInt,posAsInt+variantLength):
                        sbAtP,covAtP = parsePileupPosition(f,str(p))
                        strandBias += sbAtP
                        coverage += covAtP
                    #Normalize
                    strandBias /= variantLength
                    coverage /= variantLength
                    strandBiases.append(strandBias)
                    coverageValues.append(coverage)
        filenames = diffs[pos][variant]
        diffs[pos][variant] = {}
        if (len(strandBiases) == 0 or len(coverageValues) == 0):
            diffs[pos][variant] = 'undetermined'
            continue
        diffs[pos][variant]['strandBias'] = str(sum(strandBiases)/len(strandBiases))+':'+str(strandBiases)
        diffs[pos][variant]['coverage'] = str(sum(coverageValues)/len(coverageValues))+':'+str(coverageValues)
        diffs[pos][variant]['files'] = filenames
json.dump(diffs,open(snakemake.output[0],'w'))
