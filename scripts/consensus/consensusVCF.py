from scipy.stats import binom
import vcfpy
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure

from shared import *
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

    rcov = float(record.INFO["RCOV"])
    vcov = float(record.INFO["VCOV"])
    vcov_forward = float(record.INFO["VCOVF"])
    vcov_reverse = float(record.INFO["VCOVR"])
    tcov = rcov + vcov

    if vcov < th_cov:
        continue
    elif vcov < th_sb_cov: #use pVal
        pval = binom.pmf(vcov_forward, vcov,0.5)
        if pval > th_sb_pval:
            pass
        else:
            continue # discard the record
    else: #ratio test
        ratio = min(vcov_forward, vcov_reverse) / tcov
        if ratio > th_sbiais:
            pass
        else:
            continue # discard the record

    if vcov <= 0 or rcov + vcov <= 0:
        continue

    varRatio = vcov / tcov

    if th_het <= varRatio <= 1-th_het:
        #SNPs get the ambiguous base chars
        if len(alt) == 1 and len(ref) == 1:
            record.ALT[0].value = ambiguityLetters[frozenset({record.ALT[0].value,record.REF})]
        else:
            record.INFO['HSV'] = True

    writer.write_record(record)
