from Bio import SeqIO
from shared import *
import matplotlib.pyplot as plt
import json

x = []
y = []

uniqueKmersPerK = {}

# Read reference
reference = SeqIO.read(snakemake.input[0],'fasta')

ks = snakemake.params['krange']

for k in ks:

    # Create an empty dictionary to score kmer counts (how often each kmer is observed)
    kmers_reference = {}

    parseKmers(kmers_reference,reference,k)
    
    
    ambiguousKmers = sum(1 for kmer in kmers_reference if kmers_reference[kmer] > 1)
    totalKmers = sum(1 for kmer in kmers_reference)
    ratio = ambiguousKmers/totalKmers
    
    uniqueKmersPerK[k] = totalKmers-ambiguousKmers
    

    
    duplicates = sum(1 for kmer in kmers_reference if kmers_reference[kmer] == 2)
    
    print(k,totalKmers,duplicates)
    
    x.append(k)
    y.append(ratio)
  
# plot tobigram
  
plt.rcParams['figure.figsize'] = [10, 5]
plt.plot(x,y)
plt.xlabel('k')
plt.ylabel('fraction of ambiguous k-mers')

plt.savefig(snakemake.output['tobigram'])
with open(snakemake.output['uniqueKmers'],'w') as outfile:
	json.dump(uniqueKmersPerK,outfile)
    
