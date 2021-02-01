import sys
import altair as alt
import pandas as pd
import vcfpy
from functools import reduce

tuples = []

ambiguousLetters = ['R','Y','S','W','K','M']
iterator = iter(snakemake.input)
for pancovF,ivarF in zip(snakemake.input['pancov'],snakemake.input['ivar']):
    reader = vcfpy.Reader.from_path(pancovF)
    for record in reader:
        isHet = False
        #print(record.ALT[0].value)
        for l in ambiguousLetters:
            if l in record.ALT[0].value:
                isHet = True
                break
        if isHet:
            alleleFrequency = record.INFO['VCOV']/(record.INFO['VCOV']+record.INFO['RCOV'])
            tuples.append((pancovF,record.POS,'pancov',alleleFrequency))           
            ivartable = open(ivarF,'r').read().splitlines()[1:]
            illuminafreq = -1
            for il in ivartable:
                d = il.split()
                print(d)
                if int(d[1]) == int(record.POS):
                    illuminafreq = float(d[9])
                    tuples.append((pancovF,record.POS,'illumina',illuminaFreq))
                    break
            

    #break #remove

df = pd.DataFrame(tuples,columns=['file','pos','method','rvt'])

charts = []

for f in df['file'].unique():

    charts.append(alt.Chart(df[df.file == f],title=f).mark_rect().encode(
        y = 'pos:O',
        x = alt.X('rvt:Q',scale=alt.Scale(domain=[0,1])),
        color= 'method:N',
        tooltip = ['rvt']
    ).interactive())

concat = reduce(alt.vconcat,charts)
concat.save(snakemake.output[0])
