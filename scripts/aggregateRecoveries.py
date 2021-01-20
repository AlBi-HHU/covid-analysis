
with open(snakemake.output[0],'w') as outfile:
	total = 0
	recoveredCount = 0
	for f in snakemake.input:
		outfile.write(f+'\n')
		with open(f,'r') as infile:
			ll = infile.read().splitlines()[1:]
			for l in ll:
				outfile.write(l+'\n')
				d = l.split()
				recovered = eval(d[3])
				if recovered == -1:
					continue
				elif recovered == True:
					recoveredCount += 1
				total += 1

	outfile.write('recovered {} of {} found variants'.format(recoveredCount,total))

