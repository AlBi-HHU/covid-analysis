def parseKmers(kmers,sequence,k):
	if len(sequence) < k:
		return
	for i in range(len(sequence)+1-k):
		kmer = str(sequence[i:i+k])
		if kmer in kmers:
			kmers[kmer] += 1
		else:
			kmers[kmer] = 1
