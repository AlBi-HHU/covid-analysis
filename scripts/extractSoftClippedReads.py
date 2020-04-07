import pysam
from Bio import Seq,SeqIO,SeqRecord

alignment = pysam.AlignmentFile(snakemake.input['alignment'],'rb')
records = []

for segment in alignment.fetch():
	seq = Seq.Seq(segment.query_alignment_sequence)
	rec = SeqRecord.SeqRecord(seq, segment.qname, "", "")
	records.append(rec)

with open(snakemake.output[0],'w') as outfile:
	SeqIO.write(records,outfile,'fasta')
