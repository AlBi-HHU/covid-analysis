
from collections import Counter

import pysam

from Bio import SeqIO

def main(ref, mapping, th, consensus):

    record = next(SeqIO.parse(ref, "fasta"))
    reference = record.seq
    
    mapping = pysam.AlignmentFile(mapping, "rb")
    
    out = list()
    
    for col in mapping.pileup():
        count = Counter()
        for nuc in col.pileups:
            if nuc.is_del:
                count["D"] += 1
            elif nuc.is_refskip:
                count["I"] += 1
            else:
                count[nuc.alignment.query_sequence[nuc.query_position]] += 1

        ref_nuc = reference[col.pos]

        ref_count = count[ref_nuc]
        other_count = sum([count[k] for k in set(count.keys()) - set(ref_nuc)])

        if other_count == 0 or ref_count / other_count > th:
            out.append(ref_nuc)
        else:
            print(f"ref_nuc {ref_nuc} count {count} pos {col.pos}")
            out.append('N')

    print(f">{record.id}\n{''.join(out)}", file=open(consensus, "w"))

if "snakemake" in locals():
    main(snakemake.input["reference"], snakemake.input["mapping"], float(snakemake.params["threshold"]), snakemake.output["consensus"])
else:
    import sys

    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4])

