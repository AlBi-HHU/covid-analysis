from scipy.stats import poisson
from scipy.misc import logsumexp
import numpy as np
from math import log

#simple naive likelihood model, poisson based
def calculateLikelihoods(observations,referenceEmissions,variantEmissions,splitsteps):

    splits = np.arange(0,1,splitsteps)

    splitLikelihoods = {}
    details = {}

    for split in splits:
        #split is defined by reference frequency
        referenceFrequency = split
        variantFrequency = 1-split

        likelihood = 0

        details[split] = []

        for k in observations:
            if not k in variantEmissions:
                #not in shared set of homozygous samples, ignore?
                continue
            power = observations[k]
            probability = referenceEmissions[k]*referenceFrequency+variantEmissions[k]*variantFrequency
            #likelihood += log(probability**power)
            likelihood += power*log(probability)

        
        splitLikelihoods[split] = likelihood


    #normalize
    totalVals = np.fromiter(splitLikelihoods.values(),dtype=float)
    totalLikelihood = logsumexp(totalVals)

    for split in splitLikelihoods:
        splitLikelihoods[split] -= totalLikelihood
        #return from logspace
        splitLikelihoods[split] = np.exp(splitLikelihoods[split])#,details[split])

    return splitLikelihoods
