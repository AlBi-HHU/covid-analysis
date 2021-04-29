import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure

import seaborn as sns
import matplotlib.pyplot as plt
from shared import *

plt.rcParams["agg.path.chunksize"] = 10000

x = []
y = []

for f in snakemake.input:
    print("processing {}".format(f))
    pileup = parsePileup(f)
    for pos in pileup:
        coverage = sum(pileup[pos].values())
        block = int(pos / 100)
        x.append("{}-{}".format(block * 100, (block + 1) * 100))
        y.append(coverage)

plt.figure(figsize=(40, 12))

sns.lineplot(x, y)

plt.xlabel("position")
plt.xticks(rotation="vertical")

plt.ylabel("coverage")
plt.savefig(snakemake.output[0])
