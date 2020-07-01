import json
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import nbinom
likelihoods = json.load(open(snakemake.input['likelihoods'],'r'))
pileup = parsePileup(snakemake.input['pileup'])
reference = SeqIO.read(snakemake.input['reference'],'fasta')

def writeCall(outfile,chrom='?',pos='?',id='?',ref='?',alt='?',qual='?',filt='PASS',info='?',form='?',sample='?'):
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom,pos,id,ref,alt,qual,filt,info,form,sample))

with open(snakemake.output['vcf'],'w') as outfile:
    
    het = False

    maxLikelihoods = {}

    outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    for pos in likelihoods:
            
        #check for undetermined positions
        if not isinstance(likelihoods[pos],dict):
            continue

        for variant in likelihoods[pos]:


            #determine max

            maxSplit = '?'
            maxValue = -1.0
            #print(likelihoods[pos])
            for split in likelihoods[pos][variant]:
                val = float(likelihoods[pos][variant][split])
                if val > maxValue:
                    maxValue = val
                    maxSplit = float(split)
                

            #determine call type
            calltype = '?'

            if maxSplit < 0.2:
                writeCall(outfile,pos=pos,alt=variant,sample='1/1',form='GT',ref=reference[int(pos)-1],info='maxLikelihood={}'.format(maxSplit))
            elif maxSplit < 0.8:
                writeCall(outfile,pos=pos,alt=variant,sample='0/1',form='GT',ref=reference[int(pos)-1],info='maxLikelihood={}'.format(maxSplit))
                het = True
            else:
                pass
            maxLikelihoods[pos] = maxSplit

    fig,axes = plt.subplots(nrows=2)
    if het:
        x = []
        y = []
        bestX = 0
        bestY = float('inf')
        for split in np.arange(0,0.5,snakemake.config['splitsteps']):
            distance = 0
            for pos in maxLikelihoods:
                #find closest point
                distance += min(
                    abs(maxLikelihoods[pos]-1)**2,
                    abs(maxLikelihoods[pos]-0)**2,
                    abs(maxLikelihoods[pos]-split)**2,
                    abs(maxLikelihoods[pos]-(1-split))**2
                )
            x.append(split)
            y.append(distance)
            if distance < bestY:
                bestX = split
                bestY = distance

        axes[0].plot(x,y)
        axes[0].axvline(bestX)
        axes[0].set_xlim(left=0,right=0.5)
        axes[0].set_xlabel('split')
        axes[0].set_ylabel('squared distance over all positions')

    axes[1].hist(maxLikelihoods.values())
    axes[1].set_xlim(left=0,right=1)
    axes[1].set_xlabel('split with highest likelihood')
    axes[1].set_ylabel('count')
    fig.tight_layout(pad=3.0)
    plt.savefig(snakemake.output['mix'])
