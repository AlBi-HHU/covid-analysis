
import vcfpy
from shared import *
import sys
from copy import *


def strand_biais_filter(data, key, th_sbiais,th_sb_cov,th_sb_pval):
    total = data[False][key] + data[True][key];
    if total == 0:
        return 0

    #Decide whether to apply p-val or ratio filter
    if total < th_sb_cov: #use p-Val
        pval = binom.pmf(data[True][key],total,0.5)
        if pval > th_sb_pval:
            return total
        else:
            return 0
    else:
        ratio = min(data[False][key] / total, data[True][key] / total)
        if ratio > th_sbiais:
            return total
        else:
            return 0

degenerate = {
    frozenset(('A', 'G')): 'R',
    frozenset(('C', 'T')): 'Y',
    frozenset(('G', 'C')): 'S',
    frozenset(('A', 'T')): 'W',
    frozenset(('G', 'T')): 'K',
    frozenset(('A', 'C')): 'M',
    frozenset(('C', 'G', 'T')): 'B',
    frozenset(('A', 'G', 'T')): 'D',
    frozenset(('A', 'C', 'T')): 'H',
    frozenset(('A', 'C', 'G')): 'V',
    frozenset(('A', 'C', 'T', 'G')): 'N'
}

inv_ambiguous = {v: k for k, v in degenerate.items()}

th_sbiais = snakemake.params["th_sbiais"],
th_cov = snakemake.params["th_cov"],
th_sb_cov = snakemake.params['th_sb_cov'],
th_sb_pval = snakemake.params['th_sb_pval'],
th_het = snakemake.params["th_het"]

pileup = parsePileupStrandAwareLight(snakemake.input['pileup'])

reader = vcfpy.Reader.from_path(snakemake.input['vcf'])
writer = vcfpy.Writer.from_path(snakemake.output['vcf'], reader.header)

for record in reader:
    #We only have single variants
    ref = record.REF
    alt = record.ALT[0].value #therefore picking the first one picks the only existing variant
    pos = record.POS
    #Check for SBIAS
    upperAlt = alt.replace('-','(')
    lowerAlt = alt.lower().replace('-',')')
    if (not lowerAlt in pileup[pos]) and (not upperAlt in pileup[pos]):
        print('Can\'t find the alt allele in the pileup for either plus/minus strand: {} / {}'.format(lowerAlt,upperAlt))
        sys.exit(-1)