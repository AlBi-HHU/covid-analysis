from pysam import VariantFile
from shared import *
import json

positions = []

for labelFilePath in snakemake.input:
	with open(labelFilePath,'r') as infile:
		labels = json.load(infile)
		for position in labels:
			for allele in labels[position]:
				if labels[position][allele][0] == "HOM":
					if not position in positions:
						positions.append(position)

with open(snakemake.output[0],'w') as outfile:
	json.dump(positions,outfile)
