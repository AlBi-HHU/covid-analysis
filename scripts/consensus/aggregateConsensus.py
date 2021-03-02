from Bio import SeqIO

with open(snakemake.output[0],'w') as outfile:
	records = []

	for f in snakemake.input:
		d = f.split('/')
		id = d[-2] + ':' + d[-1]
		rec = SeqIO.read(f,'fasta')
		rec.id = id
		records.append(rec)
	SeqIO.write(records,outfile,'fasta')