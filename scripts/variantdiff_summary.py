from shared import *

iterator = iter(snakemake.input)


with open(snakemake.output[0],'w') as outfile:

    totalChanged = 0
    totalMissed = 0
    totalNew = 0

    outfile.write('{}\t{}\t{}\t{}\n'.format('POS','PANCOV','COMPARISON','PILEUP'))

    for pancovFilePath in iterator:
        comparisonFilePath = next(iterator)
        pileupFilePath = next(iterator)

        pileup = parsePileup(pileupFilePath)

        with open(pancovFilePath,'r') as pcf, open(comparisonFilePath,'r') as cof:

            comparisonData = {x.split()[0] : (x.split()[1] , x.split()[2]) for x in cof.read().splitlines()}
            pancovData = {x.split()[0] : (x.split()[1] , x.split()[2]) for x in pcf.read().splitlines()}

            for compPosition in comparisonData:

                for pancPosition in pancovData:

                    if compPosition == pancPosition:

                        pancAlt = pancovData[pancPosition][1]
                        compAlt = comparisonData[compPosition][1]

                        pancRef = pancovData[pancPosition][0]
                        compRef = comparisonData[compPosition][0] #Should be equal, add sanity check?

                        if pancAlt == compAlt:
                            # same in both vcfs
                            break
                        else:
                            # changed
                            totalChanged += 1
                            outfile.write(
                                '{}\t{}\t{}\t{}\n'.format(
                                    pancPosition,
                                    '{}->{}'.format(pancRef,pancAlt),
                                    '{}->{}'.format(compRef,compAlt),
                                    pileup[
                                        pancPosition] if pancPosition in pileup else 'no pileup available for this position'
                                )
                            )
                        break
                else:
                    totalMissed += 1
                    outfile.write(
                        '{}\t{}\t{}\t{}\n'.format(
                            compPosition,
                            'Missing',
                            '{}->{}'.format(compRef,compAlt),
                            pileup[compPosition] if compPosition in pileup else 'no pileup available for this position'
                        )
                    )
            # Check for new variants (exclusively detected by pancov)
            for pancPosition in pancovData:
                for compPosition in comparisonData:
                    if pancPosition == compPosition: #already processed
                        break
                else:
                    pancAlt = pancovData[pancPosition][1]
                    pancRef = pancovData[pancPosition][0]

                    totalNew += 1
                    outfile.write(
                        '{}\t{}\t{}\t{}\n'.format(
                            compPosition,
                            '{}->{}'.format(pancRef, pancAlt),
                            'Missing',
                            pileup[pancPosition] if pancPosition in pileup else 'no pileup available for this position'
                        )
                    )

    outfile.write(
            'New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew, totalChanged, totalMissed))
