import json
import operator
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np

pdf = matplotlib.backends.backend_pdf.PdfPages(snakemake.output['figures'])

positions = {}

for l in snakemake.input['likelihoods']:
    likelihoods = json.load(open(l,'r'))
    for pos in likelihoods:
        if not isinstance(likelihoods[pos],dict):
            continue
        if not pos in positions:
            positions[pos] = []
        for variant in likelihoods[pos]:
            maxSplit = '?'
            maxValue = -1.0
            for split in likelihoods[pos][variant]:
                val = likelihoods[pos][variant][split]
                if val > maxValue:
                    maxValue = val
                    maxSplit = split
            positions[pos].append(float(maxSplit))


for p in sorted([int(x) for x in positions.keys()]):
    p = str(p)
    fig = plt.figure()
    plt.hist(positions[p],bins=np.arange(0,1,snakemake.config['splitsteps']))
    plt.title('Position: {}'.format(p))
    plt.xlabel('Highest Likelihood')
    plt.ylabel('Count')
    pdf.savefig(fig)
    plt.close()

pdf.close()
