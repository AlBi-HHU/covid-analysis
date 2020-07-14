
from collections import Counter

import pysam

from Bio import SeqIO

degenerate = {
    frozenset(('A')): 'A',
    frozenset(('C')): 'C',
    frozenset(('G')): 'G',
    frozenset(('T')): 'T',
    frozenset(('A', 'G')): 'R',
    frozenset(('C', 'T')): 'Y',
    frozenset(('G', 'C')): 'S',
    frozenset(('A', 'T')): 'W',
    frozenset(('G', 'T')): 'K',
    frozenset(('A', 'C')): 'M ',
    frozenset(('C', 'G', 'T')): 'B',
    frozenset(('A', 'G', 'T')): 'D',
    frozenset(('A', 'C', 'T')): 'H',
    frozenset(('A', 'C', 'G')): 'V',
    frozenset(('A', 'C', 'T', 'G')): 'N'
}

def main(ref, mapping, th_cov, th_het, consensus):

    record = next(SeqIO.parse(ref, "fasta"))
    reference = record.seq
    
    mapping = pysam.AlignmentFile(mapping, "rb")
    
    out = list()
    
    for col in mapping.pileup():
        counts_in_del = Counter()
        for nuc in col.pileups:
            if nuc.is_del:
                counts_in_del["D"] += 1
            elif nuc.is_refskip:
                counts_in_del["I"] += 1
            else:
                counts_in_del[nuc.alignment.query_sequence[nuc.query_position]] += 1

        ref_nuc = reference[col.pos]
        counts = {
            'A': counts_in_del['A'],
            'C': counts_in_del['C'],
            'T': counts_in_del['T'],
            'G': counts_in_del['G']
        }

        all_count = sum(counts.values())

        if all_count < th_cov:
            out.append('N')
        elif counts[ref_nuc] / all_count > th_het:
            out.append(ref_nuc)
        else:
            print(f"ref {ref_nuc} counts {counts} all_count {all_count} {counts[ref_nuc] / all_count}")

            nucs = frozenset((k for k, v in counts.items() if v / all_count > 1 - th_het))
            print(f"degenerate {degenerate[nucs]}")
            out.append(degenerate[nucs])

    print(f">{record.id}\n{''.join(out)}", file=open(consensus, "w"))

if "snakemake" in locals():
    main(snakemake.input["reference"], snakemake.input["mapping"], int(snakemake.params["th_cov"]), float(snakemake.params["th_het"]), snakemake.output["consensus"])
else:
    import sys

    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5])

