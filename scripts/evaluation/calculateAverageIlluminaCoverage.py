import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import parsePileupStrandAwareLight


pileup = parsePileupStrandAwareLight(snakemake.input[0])
coverages = []

for pos in range(1,snakemake.config['parsePileupStrandAwareLight']):
	coverages.append(parsePileupStrandAwareLight[pos] if pos in parsePileupStrandAwareLight else 0)

with open(snakemake.output[0],'w') as outfile:
	outfile.write(sum(coverages)/len(coverages))