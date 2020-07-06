def parseKmers(kmers,sequence,k):
    if len(sequence) < k:
        return
    for i in range(len(sequence)+1-k):
        kmer = str(sequence[i:i+k])
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1


def squashStrandedness(sign):
    if sign == ')':
        return '('
    else:
        return sign.upper()

#Parse pileup analysis file to a pythonesque structure
def parsePileup(pileupfile):
    pileupAnalysis = {}
    with open(pileupfile,'r') as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            pos = int(columns[0])
            ref = columns[1]
            spread = columns[2]
            spreadDict = {}
            for val in spread.split(';'):
                entry = val.split('=')
                spreadDict[squashStrandedness(entry[0])] = int(entry[1])
            pileupAnalysis[pos] = spreadDict
    return pileupAnalysis

def parsePileupPosition(f,pos):
    with open(f,'r') as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            if pos != columns[0]:
                continue
            spread = columns[2]
            spreadDict = {}
            for val in spread.split(';'):
                entry = val.split('=')
                spreadDict[entry[0]] = int(entry[1])
            bias = 0
            for entry in spreadDict:
                if entry.islower() or entry == ')':
                    bias -= spreadDict[entry]
                else:
                    bias += spreadDict[entry]

            return bias,sum(spreadDict.values())
    
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))
