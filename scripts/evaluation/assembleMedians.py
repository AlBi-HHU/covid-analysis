import statistics
import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import *
import json

# check cohort frequencies on
meanPileups = {}

for f in snakemake.input:
    print("processing pileups: {}".format(f))
    fp = parsePileupStrandAware(f)
    for p in fp:
        if not p in meanPileups:
            meanPileups[p] = {}
        total = sum(fp[p][0].values())
        for k, v in fp[p][0].items():
            if not k in meanPileups[p]:
                meanPileups[p][k] = []
            meanPileups[p][k].append(v / total)

for p in meanPileups:
    print("calculating median values ... (pos {})".format(p))
    for k in meanPileups[p]:
        meanPileups[p][k] = statistics.median(meanPileups[p][k])

json.dump(meanPileups, open(snakemake.output[0], "w"))
