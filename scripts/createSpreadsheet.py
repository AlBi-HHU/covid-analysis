import pandas as pd
from Bio import SeqIO

mappedConsensus = pd.read_excel(snakemake.input['curation'])


letters = set()


refSeq = SeqIO.read(snakemake.input['ref'],'fasta')

consensus = {}
oldPositions = set(mappedConsensus.columns[5:-1])

changedPositions = set()
mutations = {} 
ns = {}
vs = {}
hs = {}

nrOfColumns = len(mappedConsensus.columns)-1

skip = False
for row in mappedConsensus.itertuples():
    if not skip or pd.isnull(row[1]):
        skip = True
        continue
    
    run = row[1]
    barcode = row[2]

    print('Processing dataset: {} {}'.format(run,barcode))
    
    n_count = row[3]
    v_count = row[4]
    h_count = row[5]
    variants = {}
    for columnIndex in range(6,nrOfColumns+1):
        column = mappedConsensus.columns[columnIndex-1]
        value = row[columnIndex]
        if refSeq[column-1] != value:
            variants[column] = value

        
    consensus[(run,barcode)]={
        'n_count' : n_count,
        'v_count' : v_count,
        'h_count' : h_count,
        'variants' : variants 
    }
    
    prefix = '../covid/data/auxiliary/consensus/'+mode+'/'+run+'/'+ ('0' +  str(barcode) if barcode < 10 else str(barcode))
   
    
    newSeq = SeqIO.read(prefix+'/variant_applied.fasta','fasta')
    
    mutations[(run,barcode)] = {}
    ns[(run,barcode)] = 0
    vs[(run,barcode)] = 0
    hs[(run,barcode)] = 0
    for pos in range(len(refSeq)):
        if refSeq[pos] == newSeq[pos]:
            continue
        else:
            changedPositions.add(pos+1) #add to set
            letters.add(newSeq[pos])
            if newSeq[pos] == 'N':
                ns[(run,barcode)] += 1
            else:
                if not newSeq[pos] in ['A','C','G','T']:
                    hs[(run,barcode)] += 1
            vs[(run,barcode)] += 1
            mutations[(run,barcode)][pos+1] = newSeq[pos]

            
tuples = [
        (
        'REF',
        'REF',
        'REF',
        'REF',
        'REF',
        'REF',
        'REF',
        'REF'
        )+tuple(refSeq[pos-1] for pos in allPositions)
]

def mutationEncoding(k,pos):
    ret = ''
    if pos in mutations[k]: #is mutated according to our genotyping
        ret += mutations[k][pos]   
        if pos in consensus[k]['variants']:
            if consensus[k]['variants'][pos] != mutations[k][pos]: #variant exists but differently
                ret += '(C)'
        else:
            ret += '(N)'
    else:
        ret += refSeq[pos-1]
        if pos in consensus[k]['variants']:
            ret += '(M)'
    return ret
            
for k,v in mutations.items():
    run = k[0]
    barcode = k[1]
    print('Processing:',run,barcode)
    tuples.append(
        (
        run,
        barcode,
        ns[(run,barcode)],
        consensus[(run,barcode)]['n_count'],
        vs[(run,barcode)],
        consensus[(run,barcode)]['v_count'],
        hs[(run,barcode)],
        consensus[(run,barcode)]['h_count']
        )+tuple(mutationEncoding(k,pos) for pos in allPositions)
    )  

df = pd.DataFrame(tuples, columns =['Run', 'Barcode', 'N', 'N Artic', 'Vars', 'Vars Artic', 'Hets', 'Hets Artic']+[str(pos) if pos in oldPositions else str(pos)+'(!)' for pos in allPositions]) 
df.to_excel(mode+'.xlsx')          
