import os
from shared import *
import pandas as pd
import altair as alt

print('Identifying var positions ...')

positions = set()
with open(snakemake.input['verification'],'r') as varfile:
	for l in varfile.read().splitlines():
		d = l.split('\t')
		if len(d) > 1:
			positions.add(d[0])
		else:
			continue

print('Collecting data ...')

tuples = []


for f in snakemake.input['iteratorList']:
	pileup = parsePileupStrandAwareLight(f)
	for pos in pileup:
		if pos in positions:
			tuples.append((pos,f,sum(pileup[pos].values())))

df = pd.DataFrame(tuples,columns=['pos','file','cov'])
print('plotting')

for pos in positions:
	base = alt.Chart(df[df.pos == pos])

	bar = base.mark_bar().encode(
		x=alt.X('cov:Q', bin=True),
		y='count()'
	)

	rule = base.mark_rule(color='red').encode(
		x='mean(cov):Q',
		size=alt.value(5)
	)
	os.makedirs(snakemake.output[0],exist_ok=True)
	(bar + rule).save(os.path.join(snakemake.output[0],'{}.cov.html'.format(pos)))