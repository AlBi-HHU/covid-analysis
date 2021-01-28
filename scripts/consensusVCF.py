from scipy.stats import binom
import vcfpy
from shared import *
import sys
import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)

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

th_sbiais = float(snakemake.params["th_sbiais"])
th_cov = int(snakemake.params["th_cov"])
th_sb_cov = int(snakemake.params['th_sb_cov'])
th_sb_pval = float(snakemake.params['th_sb_pval'])
th_het = float(snakemake.params["th_het"])

pileup = parsePileupStrandAwareLight(snakemake.input['pileup'])

reader = vcfpy.Reader.from_path(snakemake.input['vcf'])

#add additional header line to mark het vars
header = reader.header
header.add_info_line({"ID": "HSV", "Type": "Flag", "Number": "1",
                      "Description": "Variant might be a heterozygous SV"})

with open(snakemake.output['nMask'],'w') as outfile:
    #First of all we identify drop-out regions where our coverage is too low to make any calls, we will mask them with N letters
    for pos in range(1,snakemake.config['ref_genome_length']+1):
        if (not (pos in pileup)) or (sum(pileup[pos].values()) < th_cov):
            outfile.write('{}\t{}\n'.format(snakemake.config['ref_genome_chr'],pos))

writer = vcfpy.Writer.from_path(snakemake.output['vcf'], header )

for record in reader:
    logging.debug('Processing record: {}'.format(record))

    #We only have single variants
    ref = record.REF
    alt = record.ALT[0].value #therefore picking the first one picks the only existing variant
    pos = record.POS

    logging.debug('Corresponding pileup record: {}'.format(pileup[pos]))


    #Check for SBIAS
    upperAlt = alt.replace('-','(')
    lowerAlt = alt.lower().replace('-',')')
    if (not lowerAlt in pileup[pos]) and (not upperAlt in pileup[pos]):
        logging.debug('Can\'t find the alt allele in the pileup for either plus/minus strand: {} / {}'.format(lowerAlt,upperAlt))
        #sys.exit(-1)
    else:
        upperCount = pileup[pos][upperAlt] if upperAlt in pileup[pos] else 0
        lowerCount = pileup[pos][lowerAlt] if lowerAlt in pileup[pos] else 0
        totalCount = upperCount + lowerCount

        logging.debug('Calculated a total ALT count on both strands of: {}, comparing to threshold: {}'.format(totalCount,th_cov))

        if totalCount < th_cov:
            continue # discard the record
        elif totalCount < th_sb_cov: #use pVal
            pval = binom.pmf(upperCount,totalCount,0.5)
            if pval > th_sb_pval:
                pass
            else:
                continue # discard the record
        else: #ratio test
            ratio = min(lowerCount, upperCount) / totalCount
            if ratio > th_sbiais:
                pass
            else:
                continue # discard the record
        #Substitute heterozygosity characters on demand and write the records
        upperCount_ref = pileup[pos][ref] if ref in pileup[pos] else 0
        lowerCount_ref = pileup[pos][ref.lower()] if ref.lower() in pileup[pos] else 0
        totalCount_ref = upperCount_ref + lowerCount_ref
        varRatio = totalCount/(totalCount_ref+totalCount)

        if th_het <= varRatio <= 1-th_het:
            #SNPs get the ambiguous base chars
            if len(alt) == 1 and len(ref) == 1:
                record.ALT[0].value = degenerate[frozenset({record.ALT[0].value,record.REF})]
            else:
                record.INFO['HSV'] = True

        writer.write_record(record)