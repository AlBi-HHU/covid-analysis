from pysam import VariantFile
from shared import *
import json


outlist = []

for labelFilePath in snakemake.input:
	with open(labelFilePath,'r') as infile:

		identifier = labelFilePath #TODO: use a more elegant, easy-to-read identifier

		sampleLabel = "HOM"


		labels = json.load(infile)
		for position in labels:
			for allele in labels[position]:
				if labels[position][allele][0] == "HET":
					if sampleLabel == "HOM":
						sampleLabel = "HET"
						outlist.append((identifier,"HET"))
					#Parse frequencies

					simpleDict = {}
					strandBias = 0

					spreadDict = labels[position][allele][2]

					plus = 0
					minus = 0
					for c in spreadDict:
						if c.islower():
							minus += spreadDict[c]
						elif c.isupper():
							plus += spreadDict[c]
						if c.lower() in simpleDict:
							simpleDict[c.lower()] += spreadDict[c]
						else:
							simpleDict[c.lower()] = spreadDict[c]
					strandBias = plus-minus
					#normalize to attain frequencies
					frequencies = {k: v for k, v in sorted(simpleDict.items(), key=lambda item: item[1],reverse=True)}
					frequenciesString = ''
					for k,v in frequencies.items():
						frequenciesString += k + ':' + str(v) +' '
					outlist.append(('',str(position),str(strandBias),frequenciesString))
		if sampleLabel == "HOM":
			outlist.append((identifier,"HOM"))
						

with open(snakemake.output[0],'w') as outfile:
	for entry in outlist:
		outfile.write('\t'.join(entry)+'\n')
