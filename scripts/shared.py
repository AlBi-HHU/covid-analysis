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

#Parse pileup analysis file to a pythonesque structure
def parsePileupStrandAware(pileupfile):
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
                spreadDict[entry[0]] = int(entry[1])
            biases = []
            processed = set()
            for entry in spreadDict:
                mode = 'low'
                if entry == ')' or entry.islower():
                    mode = 'low'
                else:
                    mode = 'high'

                if entry.upper() in processed:
                    continue
                else:
                    processed.add(entry.upper())

                for other in spreadDict:
                    if other == entry:
                        continue
                    if mode == 'low' and other.lower() == entry:
                        biases += [abs(spreadDict[entry]-spreadDict[other])]
                        break
                    elif mode == 'high' and other == entry.lower():
                        biases += [abs(spreadDict[entry]-spreadDict[other])]
                        break
            pileupAnalysis[pos] = (spreadDict,biases)
    return pileupAnalysis


def parsePileupStrandAwareLight(pileupfile):
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
                spreadDict[entry[0]] = int(entry[1])
            pileupAnalysis[pos] = spreadDict
    return pileupAnalysis


#returns the strand bias and the coverage for a given position and a pileup file
def parsePileupPosition(f,pos):
    with open(f,'r') as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            if pos != columns[0]: #column 0 contains the position
                continue #we loop until we hit the position we are interested in
            spread = columns[2] #column 2 contains the observed alleles
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
        raise Exception('Can\'t find an entry for position {} in file {}'.format(pos,f))
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))

def getStrandBias(pileupForPosition,altallele):

    altallele = '(' if altallele == '-' else altallele
    plus = 0
    total = 0
    for allele,count in pileupForPosition.items():
        if (allele.isupper() or allele == '(') and allele == altallele:
            plus += count
        if squashStrandedness(allele) == altallele:
            total += count
    return plus/total if total != 0 else 0

def getCoverage(pileupForPosition,altallele):
    altallele = '(' if altallele == '-' else altallele
    total = 0
    for allele,count in pileupForPosition.items():
        if squashStrandedness(allele) == altallele:
            total += count
    return total

def getMinorStrandAbs(pileupForPosition,altallele):
    altallele = '(' if altallele == '-' else altallele
    total = 0
    for allele,count in pileupForPosition.items():
        if squashStrandedness(allele) == altallele and (allele == ')' or allele.islower()):
            total += count
    return total

def getMinorStrandFrequency(pileupForPosition,altallele):
    sb = getStrandBias(pileupForPosition,altallele)
    return min(1-sb,sb)