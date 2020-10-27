generatedDelKmers = set()

with open('nCoV-2019.reference.k11.csv','r') as infile:

    lines = list(infile.read().splitlines())

    for l in lines:
        kmer = l.split()[0]
        
        lastChar = 'N'
        curLength = 0
        startingIndex = 0
        
        for idx,c in enumerate(kmer):
            if c == lastChar:
                if curLength == 0:
                    startingIndex = idx
                curLength += 1
                if curLength >= 4:
                    delKmer = kmer[:startingIndex]+kmer[startingIndex+1:]
                    #print(kmer,delKmer)
                    generatedDelKmers.add(delKmer)
            else:
                curLength = 0
            lastChar = c
            
    for l in lines:
        kmer = l.split()[0]
        
        for delKmer in generatedDelKmers:
            if delKmer in kmer:
                #print(delKmer,kmer)
                break
        else:
            print(kmer)
