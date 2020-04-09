from pysam import VariantFile
from shared import *
import json

emissionVectors = {}
positions = []
output = {}

for labelFilePath in snakemake.input:
	with open(labelFilePath,'r') as infile:
		labels = json.load(infile)
		for position in labels:
			for allele in labels[position]:
				if labels[position][allele][0] == "HOM":
					if not position in positions:
						positions.append(position)
						emissionVectors[position] = {}
					if not allele in emissionVectors[position]:
						emissionVectors[position][allele] = []
					emissionVectors[position][allele].append(labels[position][allele][2])

for position in positions:
	output[position] = {}
	for allele in emissionVectors[position]:

		output[position][allele] = {}

		alleleSet = set.union(*[set(e.keys()) for e in emissionVectors[position][allele]])
		outdict = {
			k : 0 
			for k in alleleSet
		}
		numberOfEntries = len(emissionVectors[position][allele])
		for entry in emissionVectors[position][allele]:
			sumOfCounts = sum(entry.values())
			for k,v in entry.items():
				outdict[k] += v/sumOfCounts
		for k,v in outdict.items():
			outdict[k] = v/numberOfEntries

		
		output[position][allele] = outdict
	

with open(snakemake.output[0],'w') as outfile:
	json.dump(output,outfile)
