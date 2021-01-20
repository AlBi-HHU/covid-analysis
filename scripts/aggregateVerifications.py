
with open(snakemake.output[0],'w') as outfile:
	total = 0
	verified = 0
	for f in snakemake.input:
		outfile.write(f+'\n')
		with open(f,'r') as infile:
			ll = infile.read().splitlines()
			for l in ll:
				outfile.write(l+'\n')
				d = l.split()
				reject = d[3]
				if not reject:
					verified += 1
				total += 1

	outfile.write('{},{}'.format(total,verified))

