from shared import *
import vcfpy


def altToText(ALT):
    return '/'.join([x.value for x in ALT])

def altEqualVar(alt,var):
    if len(alt) != 1:
        return False
    else:
        if alt[0].value == var:
            return False
    return True

with open(snakemake.output[0],'w') as outfile:

    totalNew = 0
    totalChanged = 0
    totalMissed = 0

    outfile.write('{}\t{}\t{}\t{}\n'.format('Pos', 'Pancov', 'ComparedMethod','Pileup (First Pos)'))

    iterator = iter(snakemake.params['vcfs'])
    for comparison_table in iterator:

        pancov_vcf = next(iterator)
        pileup = parsePileup(next(iterator))

        if comparison_table == '':
            #TODO: Log Debug Output? ...
            print('No manual curation (GISAID Seq) found for: {}'.format(oancov_vcf))
            continue


        outfile.write(' \t{}\t{}\n'.format(pancov_vcf,comparison_table))

        originalVCF = [r for r in vcfpy.Reader(open(pancov_vcf, 'r'))] #pancov
        comparisonTable = {}
        with open(comparison_table,'r') as infile:
            for l in infile.read().splitlines():
                data = l.split()
                pos1based = int(data[0])
                ref = data[1]
                alt = data[2]
                comparisonTable[pos1based] = (ref,alt)

        processedPositions = []

        #Check all records in the vcf we compare ourselves to
        for onebasedcomp in comparisonTable:
            for originalRecord in originalVCF:
                onebasedorig = int(originalRecord.POS)#-1
                if record.POS == originalRecord.POS:
                    if altEqualVar(originalRecord.ALT,comparisonTable[onebasedcomp][1]):
                        #same in both vcfs
                        break
                    else:
                        #changed
                        totalChanged += 1
                        outfile.write(
                            '{}\t{}\t{}\t{}\n'.format(
                                record.POS,
                                '{}->{}'.format(originalRecord.REF,altToText(originalRecord.ALT)),
                                '{}->{}'.format(comparisonTable[onebasedcomp][0],comparisonTable[onebasedcomp][1]),
                                pileup[onebasedorig] if onebasedorig in pileup else 'no pileup available for this position'
                            )
                        )
                    break
            else:
                totalMissed += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        record.POS,
                        'Missing',
                        '{}->{}'.format(comparisonTable[onebasedcomp][0],comparisonTable[onebasedcomp][1]),
                        pileup[onebasedcomp] if onebasedcomp in pileup else 'no pileup available for this position'
                    )
                )
        #Check for new variants (exclusively detected by pancov)
        for originalRecord in originalVCF:
            onebasedorig = int(originalRecord.POS) #- 1
            for onebasedcomp in comparisonTable:
                if onebasedcomp == originalRecord.POS:
                    break
            else:
                totalNew += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        originalRecord.POS,
                        '{}->{}'.format(originalRecord.REF,altToText(originalRecord.ALT)),
                        'Missing',
                        pileup[onebasedorig] if onebasedorig in pileup else 'no pileup available for this position'
                    )
                )

    outfile.write('New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew,totalChanged,totalMissed))