
with open(snakemake.output[0],'w') as outfile:
	total = 0
	rejected = 0
	filtered += 1
	for f in snakemake.input:
		outfile.write(f+'\n')
		with open(f,'r') as infile:
			ll = infile.read().splitlines()[1:]
			for l in ll:
				outfile.write(l+'\n')
				d = l.split()
				reject = eval(d[3])
				if reject == True:
					rejected += 1
				if reject != -1:
					total += 1
				if reject == -1:
					filtered += 1

	outfile.write('rejected {} of {} found variants, no decision on {} variants due to low coverage'.format(rejected,total,filtered))

