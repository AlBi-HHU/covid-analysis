from shared import *

iterator = iter(snakemake.input)

for pancovFilePath in iterator:
    comparisonFilePath = iterator.next()
    pileupFilePath = iterator.next()

    pileup = parsePileup(pileupFilePath)

    with open(pancovFilePath,'r') as pcf, open(comparisonFilePath,'r') as cof,snakemake.output[0] as outfile:

        comparisonData = {x.split()[0] : (x.split()[1] , x.split()[2]) for x in cof.read().splitlines()}

        for record in newVCF:
            onebasedcomp = int(record.POS)  # -1

            # filter based on strand bias and coverage

            for originalRecord in originalVCF:

                onebasedorig = int(originalRecord.POS)  # -1

                if record.POS == originalRecord.POS:
                    if altEqual(record.ALT, originalRecord.ALT):
                        # same in both vcfs
                        break
                    else:
                        # changed
                        totalChanged += 1
                        outfile.write(
                            '{}\t{}\t{}\t{}\n'.format(
                                record.POS,
                                '{}->{}'.format(originalRecord.REF, altToText(originalRecord.ALT)),
                                '{}->{}'.format(record.REF, altToText(record.ALT)),
                                pileup[
                                    onebasedorig] if onebasedorig in pileup else 'no pileup available for this position'
                            )
                        )
                    break
            else:
                totalMissed += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        record.POS,
                        'Missing',
                        '{}->{}'.format(record.REF, altToText(record.ALT)),
                        pileup[onebasedcomp] if onebasedcomp in pileup else 'no pileup available for this position'
                    )
                )
        # Check for new variants (exclusively detected by pancov)
        for originalRecord in originalVCF:
            onebasedorig = int(originalRecord.POS)  # - 1
            for record in newVCF:
                if record.POS == originalRecord.POS:
                    break
            else:
                totalNew += 1
                outfile.write(
                    '{}\t{}\t{}\t{}\n'.format(
                        originalRecord.POS,
                        '{}->{}'.format(originalRecord.REF, altToText(originalRecord.ALT)),
                        'Missing',
                        pileup[onebasedorig] if onebasedorig in pileup else 'no pileup available for this position'
                    )
                )

    outfile.write(
        'New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew, totalChanged, totalMissed))
