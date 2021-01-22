import sys
import altair as alt
import pandas as pd

tuples = []

with open(snakemake.input[0],'r') as infile:
    sample = '???'
    for l in infile.read().splitlines():
        split = l.split('\t')
        if len(split) > 1:
            if len(split) == 7:
                tuples.append((sample,int(split[0]),split[1],split[2],split[3],split[4],split[5],split[6]))
            elif len(split) == 6:
                tuples.append((sample,int(split[0]),split[1],split[2],split[3],'',split[4],split[5]))
        else:
            sid = l.split('/')[-1].split('_')
            if len(sid) == 3:
                #print(sid)
                sample =sid[0]+'/'+sid[1]
            
df = pd.DataFrame(tuples,columns=['sample','pos','ref','alt','recovered','comment','pileupillu','pileupnano'])
#print(df['recovered'])
df['recovered'] = df['recovered'].replace({'-1' : 'AlexFiltered'})
#print(df['recovered'])


make = pd.DataFrame({'sample': list(df['sample'].unique())})
selection = alt.selection_multi(fields=['sample'])
color = alt.condition(selection, alt.value('green'), alt.value('lightgray'))
make_selector = alt.Chart(make).mark_rect().encode(x='sample',color=color).add_selection(selection)

chart = alt.Chart(df).mark_rect().encode(
    y = 'pos:O',
    x = 'sample:N',
    color= alt.Color('recovered',scale=alt.Scale(
        domain = ['True','False','AlexFiltered'],
        range=['blue','orange','grey']
    )),
    tooltip = ['ref','alt','comment','pileupillu','pileupnano']
).transform_filter(
    selection
).interactive()

alt.vconcat(make_selector,chart,padding=64).save(snakemake.output[0])
