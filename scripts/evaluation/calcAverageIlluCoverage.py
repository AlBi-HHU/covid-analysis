import sys

sys.path.append(
    "scripts"
)  # Hackfix but results in a more readable scripts folder structure
from shared import *
import pandas as pd

tuples = []

for f in snakemake.input:
    pileup = parsePileupStrandAwareLight(f)
    for pos, emissions in pileup.items():
        tuples.append((f, int(pos), sum(emissions.values())))

df = pd.DataFrame(tuples, columns=["file", "pos", "cov"]).groupby("pos").mean()
df.to_csv(snakemake.output[0])
