generatedDelKmers = set()
actualKmers = set()

with open(snakemake.input[0], "r") as infile, open(snakemake.output[0], "w") as outfile:

    # read the lines and transform to list so it can be iterated twice
    lines = list(infile.read().splitlines())

    for l in lines:
        kmer = l.split()[0]

        actualKmers.add(kmer)

        lastChar = "N"
        curLength = 0
        startingIndex = 0

        for idx, c in enumerate(kmer):
            if c == lastChar:
                if curLength == 0:
                    startingIndex = idx
                curLength += 1
                # check all homopolymers and create the kmers that result from deletion
                if curLength >= 4:
                    delKmer = kmer[:startingIndex] + kmer[startingIndex + 1 :]
                    # print(kmer,delKmer)
                    generatedDelKmers.add(delKmer)
            else:
                curLength = 0
            lastChar = c

    for delKmer in generatedDelKmers:
        for kmer in actualKmers:
            if delKmer in kmer:
                # print(delKmer,kmer)
                break
        else:
            outfile.write("{}\n".format(delKmer))
