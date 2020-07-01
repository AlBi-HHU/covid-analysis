import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import json

kmer_counts = json.load(open(snakemake.input['kmerProfile'],'r'))
uniqueKmersPerK = json.load(open(snakemake.input['uniqueKmers'],'r'))
k= snakemake.params['k']

plt.rcParams['figure.figsize'] = [16, 24]

fig,axs = plt.subplots(2)

kmer_counts_list = [kmer_counts[kmer] for kmer in  kmer_counts]
b = np.amax(kmer_counts_list)

hist,_,_ = axs[0].hist(kmer_counts_list, bins=b)
axs[0].set_title('raw k-mer histogram')

local_mins = scipy.signal.find_peaks_cwt(-hist,np.arange(30,70))
threshold = local_mins[0]


# Create a list of k-mers that occur more frequent than the cutoff-point
filteredKmers = [kmer for kmer in kmer_counts if kmer_counts[kmer] > threshold]
kmer_counts_list_filtered = [kmer_counts[kmer] for kmer in filteredKmers]
b = np.arange(np.amax(kmer_counts_list_filtered) + 1)
hist_filtered, _, _ = axs[1].hist(kmer_counts_list_filtered, bins=b)

# Compare to expected value
axs[1].set_title(
    'error-threshold: {} | k: {} | unique k-mers: {}/{} expected'.format(threshold,k,len(filteredKmers),uniqueKmersPerK[str(k)])
    ,fontsize=9
)


    
fig.text(0.5,0.04,'k-mer frequency',ha='center')
fig.text(0.04,0.5,'k-mer count',va='center',rotation='vertical')
plt.savefig(snakemake.output[0])

