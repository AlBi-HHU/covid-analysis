from scipy.stats import binom
import vcfpy
from ./../shared import *
import logging

#Enable logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)

#Reassign parameters from snakemake config for better readability and type casting (?)
th_sbiais = float(snakemake.config["consensusStrandBiais"])
th_cov = int(snakemake.config["consensusMinCov"])
th_sb_cov = int(snakemake.config['consensusPValCoverage'])
th_sb_pval = float(snakemake.config['consensusPValCutoff'])
th_het = float(snakemake.config["thresholdHomCall"])

#See shared.py
pileup = parsePileupStrandAwareLight(snakemake.input['pileup'])

reader = vcfpy.Reader.from_path(snakemake.input['vcf'])

#add additional header line to mark het vars
header = reader.header #Copy the existing header
header.add_info_line({"ID": "HSV", "Type": "Flag", "Number": "1","Description": "Variant might be a heterozygous SV"})

# First of all we identify drop-out regions where our coverage is too low to make any calls, we will mask them with N letters
with open(snakemake.output['nMask'],'w') as outfile:
    for pos in range(1,snakemake.config['ref_genome_length']+1): #Note, that coordinates are 1-based
        #Positions that are not in the pileup are not covered at all, for all others: Check the coverage and compare to the threshold
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
                record.ALT[0].value = ambiguityLetters[frozenset({record.ALT[0].value,record.REF})]
            else:
                record.INFO['HSV'] = True

        writer.write_record(record)