
from collections import Counter,defaultdict

import json

import pysam

from Bio import SeqIO



def main(ref, mapping, th_cov, th_het, th_sbiais,th_sb_cov,th_sb_pval, consensus,variant_index):

    vi = {}

    record = next(SeqIO.parse(ref, "fasta"))
    reference = record.seq
    
    mapping = pysam.AlignmentFile(mapping, "rb")

    out = ['N'] * (mapping.lengths[0] + 1)
    
    for col in mapping.pileup():
        counts_in_del = defaultdict(Counter)
        for nuc in col.pileups:
            if nuc.is_del:
                counts_in_del[None]["D"] += 1
            elif nuc.is_refskip:
                counts_in_del[None]["I"] += 1
            else:
                counts_in_del[nuc.alignment.is_reverse][nuc.alignment.query_sequence[nuc.query_position]] += 1

        ref_nuc = reference[col.pos]
        counts = {
            'A': strand_biais_filter(counts_in_del, 'A', th_sbiais,th_sb_cov,th_sb_pval),
            'C': strand_biais_filter(counts_in_del, 'C', th_sbiais,th_sb_cov,th_sb_pval),
            'T': strand_biais_filter(counts_in_del, 'T', th_sbiais,th_sb_cov,th_sb_pval),
            'G': strand_biais_filter(counts_in_del, 'G', th_sbiais,th_sb_cov,th_sb_pval),
        }

        all_count = sum(counts.values())

        if all_count < th_cov:
            #print(f"ref {ref_nuc} counts {counts} all_count {all_count}")
            out[col.pos] = 'N'
        elif counts[ref_nuc] / all_count > th_het:
            out[col.pos] = ref_nuc
        else:
            nucs = frozenset((k for k, v in counts.items() if v / all_count > 1 - th_het))

            out[col.pos] = degenerate[nucs]
            vi[col.pos] = degenerate[nucs]

            #print(col.pos, degenerate[nucs])
            #print("ref_nuc", ref_nuc)
            #print("consensus", degenerate[nucs])
            #print("all_counts", counts_in_del)
            #print("count", counts)
            

    print(f">{record.id}\n{''.join(out)}", file=open(consensus, "w"))

    json.dump(vi,open(variant_index,'w'))


    

if "snakemake" in locals():
    main(snakemake.input["reference"], snakemake.input["mapping"], int(snakemake.params["th_cov"]), float(snakemake.params["th_het"]), float(snakemake.params["th_sbiais"]),float(snakemake.params["th_sb_cov"]),float(snakemake.params["th_sb_pval"]), snakemake.output["consensus"],snakemake.output['variant_index'])
else:
    import sys

    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7]), sys.argv[8], sys.argv[9])

