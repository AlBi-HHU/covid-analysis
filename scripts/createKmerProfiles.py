from Bio import SeqIO
from shared import *
import json

scReads = SeqIO.parse(snakemake.input[0],'fasta')

k = snakemake.params['k']

kmers = {}

for read in scReads:
	parseKmers(kmers,read.seq,k)

with open(snakemake.output[0],'w') as outfile:
	json.dump(kmers,outfile)
