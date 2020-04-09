from pysam import VariantFile
from shared import *
import json

positions = []
emissionProbabilities = {}

for labelFilePath in snakemake.input:
	with open(labelFilePath,'r') as infile:
		labels = json.load(infile)
		for position in labels:
			for allele in labels[position]:
				if labels[position][allele][0] == "HOM":
					if not position in positions:
						positions.append(position)
						emissionProbabilities[position] = []
					emissionProbabilities[position].append(labels[position][allele][1])

output = {
	position : (sum(emissionProbabilities[position])/len(emissionProbabilities[position]),emissionProbabilities[position])
	for position in positions
}


with open(snakemake.output[0],'w') as outfile:
	json.dump(output,outfile)
