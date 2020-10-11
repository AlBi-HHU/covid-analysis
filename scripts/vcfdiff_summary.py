from shared import *

import vcfpy

'''
#Part 1: We want to map the files correctly
mapping = {}
for vcf in snakemake.input['vcfs']:
    mapping[vcf] = 'unknown'
    #key = pancov vcf, value = comparison vcf

for vcf in snakemake.input['comparisonVcfs']:
    parsedData = vcf.split('/')
    new_run = parsedData[-2]
    new_barcode = parsedData[-1].split('.')[0][7:]

    for originalVcf in mapping:
        parsedData = originalVcf.split('/')
        orig_run = parsedData[-3]
        orig_barcode = parsedData[-2]

        if orig_run == new_run and orig_barcode == new_barcode:
            mapping[originalVcf] = vcf
            print('mapped {} to {}'.format(vcf,originalVcf))
            break
    else:
        print('Could not map anything to: {}'.format(vcf))

for vcf in mapping:
    if mapping[vcf] == 'unknown':
        print('Could not map anything to: {} '.format(vcf))

print(mapping)
'''
#Step 2 Analysis

def altToText(ALT):
    return '/'.join([x.value for x in ALT])

def altEqual(a1,a2):
    if len(a1) != len(a2):
        return False
    else:
        for x,y in zip(a1,a2):
            if x.value != y.value:
                return False
    return True

with open(snakemake.output[0],'w') as outfile:

    totalNew = 0
    totalChanged = 0
    totalMissed = 0

    outfile.write('{}\t{}\t{}\t{}\n'.format('Pos', 'Pancov', 'ComparedMethod','Pileup (First Pos)'))

    iterator = iter(snakemake.input['vcfs'])
    for comparison_vcf in iterator:
        pancov_vcf = next(iterator)
        pileup = parsePileup(next(iterator))

        outfile.write(' \t{}\t{}\n'.format(pancov_vcf,comparison_vcf))

        originalVCF = [r for r in vcfpy.Reader(open(pancov_vcf, 'r'))] #pancov
        newVCF = [r for r in vcfpy.Reader(open(comparison_vcf, 'r'))]


        processedPositions = []

        #Check all records in the vcf we compare ourselves to
        for record in newVCF:
            for originalRecord in originalVCF:
                if record.POS == originalRecord.POS:
                    if altEqual(record.ALT,originalRecord.ALT):
                        #same in both vcfs
                        break
                    else:
                        #changed
                        totalChanged += 1
                        outfile.write(
                            '{}\t{}\t{}\t{}\n'.format(
                                record.POS,
                                '{}->{}'.format(record.REF,altToText(originalRecord.ALT)),
                                '{}->{}'.format(record.REF,altToText(record.ALT)),
                                pileup[int(originalRecord.POS)] if int(originalRecord.POS) in pileup else 'no pileup available for this position'
                            )
                        )
                    break
            else:
                totalMissed += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        record.POS,
                        'Missing',
                        '{}->{}'.format(record.REF,altToText(record.ALT)),
                        pileup[int(originalRecord.POS)] if int(originalRecord.POS) in pileup else 'no pileup available for this position'
                    )
                )

        for originalRecord in originalVCF:
            for record in newVCF:
                if record.POS == originalRecord.POS:
                    break
            else:
                totalNew += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        originalRecord.POS,
                        '{}->{}'.format(record.REF,altToText(originalRecord.ALT)),
                        'Missing',
                        pileup[int(originalRecord.POS)] if int(originalRecord.POS) in pileup else 'no pileup available for this position'
                    )
                )

    outfile.write('New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew,totalChanged,totalMissed))