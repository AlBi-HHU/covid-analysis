import pysam
from Bio import Seq, SeqIO, SeqRecord

from shared import rev_comp

alignment = pysam.AlignmentFile(snakemake.input["alignment"], "rb")
records = []

for segment in alignment.fetch():
    # Check for empty segments (Not sure why they would occur?)
    if segment.query_alignment_sequence == None:
        continue
    # Reverse if required
    if segment.is_reverse:
        seq = Seq.Seq(rev_comp(segment.query_alignment_sequence))
    else:
        seq = Seq.Seq(segment.query_alignment_sequence)

    rec = SeqRecord.SeqRecord(seq, segment.qname, "", "")
    records.append(rec)

with open(snakemake.output[0], "w") as outfile:
    SeqIO.write(records, outfile, "fasta")
