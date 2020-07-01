import json
from Bio import SeqIO


likelihoods = json.load(open(snakemake.input['likelihoods'],'r'))
