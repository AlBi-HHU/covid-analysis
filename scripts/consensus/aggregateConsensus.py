from Bio import SeqIO

with open(snakemake.output[0],'w') as outfile:
	records = []

	for f in snakemake.input:
		rec = SeqIO.read(f,'fasta')
		rec.id = f
		records.append(rec)
	SeqIO.write(records,outfile,'fasta')