import sys
import altair as alt
import pandas as pd
import vcfpy

tuples = []

ambiguousLetters = ['R','Y','S','W','K','M']

for f in snakemake.input:
    reader = vcfpy.Reader.from_path(f)
    for record in reader:
        isHet = False
        for l in ambiguousLetters:
            if l in record.ALT:
                isHet = True
                break
        if isHet:
            tuples.append((f,record.POS,record.INFO['RVT']))

df = pd.DataFrame(tuples,columns=['file','pos','rvt'])

chart = alt.Chart(df).mark_rect().encode(
    y = 'pos:O',
    x = 'file:N',
    color= 'rvt:Q'
).interactive().save(snakemake.output[0])
