
with open(snakemake.output[0],'w') as outfile:
	total = 0
	verified = 0
	for f in snakemake.input:
		with open(f,'r') as infile:
			ll = infile.read().splitlines()
			for l in ll:
				d = l.split()
				reject = d[3]
				if not reject:
					verified += 1
				total += 1

	outfile.write(total,verified)

