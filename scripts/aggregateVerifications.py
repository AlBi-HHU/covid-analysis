
with open(snakemake.output[0],'w') as outfile:
	total = 0
	rejected = 0
	for f in snakemake.input:
		outfile.write(f+'\n')
		with open(f,'r') as infile:
			ll = infile.read().splitlines()
			for l in ll:
				outfile.write(l+'\n')
				d = l.split()
				reject = eval(d[3])
				if reject:
					rejected += 1
				total += 1

	outfile.write('rejected {} of {} found variants'.format(rejected,total))

