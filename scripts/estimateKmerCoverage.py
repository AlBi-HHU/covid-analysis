from Bio import SeqIO

from itertools import chain

import numpy as np
import matplotlib.pyplot as plt

import scipy.signal

import logging

from scipy.stats import poisson
from shared import *

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format="%(asctime)s:%(levelname)s:%(message)s")


#########################################################################################

# Assign parameter k to a local variable for better readability
k = int(snakemake.params['k'])

# Create an empty dictionary to score kmer counts (how often each kmer is observed)
kmer_counts = {}

# Parse reads
reads1 = SeqIO.parse(snakemake.input['read1'], 'fastq')
reads2 = SeqIO.parse(snakemake.input['read2'], 'fastq')
all_reads = chain(reads1, reads2)

# We could use a probabilistic method (sketching) here to improve memory/runtime
# For now: Simply process all reads and identify k-mers (this requires a fair bit of memory and time)
for read in all_reads:
    parseKmers(kmer_counts, read.seq, k)  # see shared.py

# Transform into a list of counts
kmer_counts_list = [kmer_counts[kmer] for kmer in kmer_counts]

# Determine the highest value that occurs
b = np.amax(kmer_counts_list)

# Create a histogram
hist,_ = np.histogram(kmer_counts_list, bins=b)

# We assume a roughly bimodal distribution of actual k-mers and sequencing error k-mers
# Choosing local min here, see DOI: 10.1093/bib/bbv029 to create a cutoff

local_mins = scipy.signal.argrelextrema(hist, np.less)[0]
threshold = local_mins[0]

logging.info('Using threshold: {} as a cutoff-point for error-k-mers'.format(threshold))

# Create a list of k-mers that occur more frequent than the cutoff-point
kmer_counts_list_filtered = [kmer_counts[kmer] for kmer in kmer_counts if kmer_counts[kmer] > threshold]

# Calculate kmer errors
kmer_errors = (sum(kmer_counts_list) - sum(kmer_counts_list_filtered)) / sum(kmer_counts_list)

with open(snakemake.output['kmererror'], 'w') as outfile:
    outfile.write('Kmer Error:' + '\t' + str(kmer_errors) + '\n')

b = np.arange(np.amax(kmer_counts_list_filtered) + 1)

hist_filtered, _, _ = plt.hist(kmer_counts_list_filtered, bins=b)

plt.savefig(snakemake.output['histogramRaw'])
plt.clf()

poisson_lambda = np.argmax(hist_filtered)
mean = np.mean(kmer_counts_list_filtered)

with open(snakemake.output['mean'], 'w') as outfile:
    outfile.write('Mean:' + '\t' + str(mean) + '\n' + 'Lambda:' + '\t' + str(poisson_lambda))

kmer_counts_list_filtered_cutoff_large_vals = [x for x in kmer_counts_list_filtered if x < CUTOFF_FACTOR_MEAN * mean]

b = np.arange(np.amax(kmer_counts_list_filtered_cutoff_large_vals) + 1)

hist_final, _, _ = plt.hist(kmer_counts_list_filtered_cutoff_large_vals, bins=b)


# Calculate Peak
actualMax = poisson.pmf(np.argmax(hist_final),poisson_lambda)
desiredMax = np.amax(hist_final)
factor = desiredMax/actualMax

plt.plot(b, poisson.pmf(b, poisson_lambda) * factor,color='orange')
plt.axvline(mean, color='r')

plt.savefig(snakemake.output['histogram'])
