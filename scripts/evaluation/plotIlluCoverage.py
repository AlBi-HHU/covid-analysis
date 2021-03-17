import sys
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import *
import altair as alt
import pandas as pd

tuples = []

for f in snakemake.input:
    pileup = parsePileupStrandAwareLight(f)
    for pos,emissions in pileup.items():
        tuples.append(
            (f,int(pos),sum(emissions.values()))
        )

df = pd.DataFrame(tuples,columns=['file','pos','cov'])

chart = alt.Chart(df).mark_rect().encode(
    y = 'mean(cov):Q',
    x = 'pos:O',
    tooltip = ['cov','pos']
).interactive()

chart.save(snakemake.output[0])
