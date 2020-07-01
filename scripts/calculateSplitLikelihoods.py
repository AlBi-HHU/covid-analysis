import json
from probabilistic import *
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from shared import *

pdf = matplotlib.backends.backend_pdf.PdfPages(snakemake.output['figures'])


pileup = parsePileup(snakemake.input['pileup'])
emissions = json.load(open(snakemake.input['emissions'],'r'))


likelihoods = {}
for position in sorted(emissions.keys(),key= lambda x : int(x)):
    likelihoods[position] = {}
    for variant in emissions[position]:
        
        if variant == 'reference':
            continue

        idx = int(position)
        
        if not idx in pileup:
            print('Dropout at position: {}'.format(idx))
            continue

        observedEmissions = pileup[idx]
        '''
        if not position in emissions:
            likelihoods[position] = '??? (No HOM/reference emissions known at given position)'
            continue
            #raise Exception('No reference emissions recorded for position: {}'.format(position))
        '''

        likelihoods[position][variant] = calculateLikelihoods(observedEmissions,emissions[position]['reference'],emissions[position][variant],snakemake.config['splitsteps'])

        data = sorted(likelihoods[position][variant].items())
        x,y = zip(*data)
        fig = plt.figure()
        plt.plot(x,y)
        plt.ylim(0,1)
        plt.xlim(0,1)
        plt.title('Position: {} Variant: {}'.format(position,variant))
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        pdf.savefig(fig)
        plt.close()

with open(snakemake.output['likelihoods'],'w') as outfile:
    json.dump(likelihoods,outfile)

pdf.close()
