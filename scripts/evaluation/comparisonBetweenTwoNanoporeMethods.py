def processData(file):

	#a dictionary will be the return value
	ret = {}

	#split the file based on linebreaks
	lines = file.read().splitlines()
	currentPosition = 0
	while (True):
		#Fetch a new section that relates to one sample
		data = lines[currentPosition].split('\t')
		if data[0].startswith("Real"): #We reached the final stat section of the file
			break
		else:
			sample = data[0]
			ret[sample] = {}
			currentPosition += 2 #skip the header section
			while (True):
				data = lines[currentPosition].split('\t')
				#escape condition: we reached another sample
				if len(data) == 1: #This is another header for the next sample
					break
				else:
					ret[sample][data[0]] = data[3] #key: position, value: nanopore variant call
					currentPosition += 1
	return ret



with open(snakemake.input['method1'],'r') as meth1file, \
	 open(snakemake.input['method2'],'r') as meth2file, \
	 open(snakemake.output[0],'w') as outfile:

	meth1data = processData(meth1file)
	meth2data = processData(meth2file)
	for sample in set(meth1data.keys()).union(meth2data.keys()):
		outfile.write('Comparison: {} \n'.format(sample))
		outfile.write('{}\t{}\t{}\n'.format('Position',snakemake.params['method1name'],snakemake.params['method2name']))
		if sample not in meth1data:
			outfile.write('Not found in {} \n'.format(snakemake.params['method1name']))
			for k,v in meth2data[sample]:
				outfile.write('{}\t{}\n'.format(k,v))
		if sample not in meth2data:
			outfile.write('Not found in {} \n'.format(snakemake.params['method2name']))
			for k, v in meth1data[sample]:
				outfile.write('{}\t{}\n'.format(k, v))
		for k in set(meth1data[sample].keys()).union(set(meth2data[sample].keys())):
			if k not in meth1data[sample]:
				outfile.write('{}\t{}\t{}\n'.format(k,'Not called',meth2data[sample][k]))
			elif k not in meth2data[sample]:
				outfile.write('{}\t{}\t{}\n'.format(k, meth1data[sample][k], 'Not called'))
			elif meth1data[sample][k] != meth2data[sample][k]:
				outfile.write('{}\t{}\t{}\n'.format(k, meth1data[sample][k], meth2data[sample][k]))
