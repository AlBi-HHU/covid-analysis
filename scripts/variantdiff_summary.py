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

        pileup = parsePileupStrandAware(pileupFilePath)

        with open(pancovFilePath,'r') as pcf, open(comparisonFilePath,'r') as cof:

            outfile.write('{}\t{}\n'.format(pancovFilePath,comparisonFilePath))

            comparisonData = {x.split()[0] : (x.split()[1] , x.split()[2]) for x in cof.read().splitlines()}
            pancovData = {x.split()[0] : (x.split()[1] , x.split()[2]) for x in pcf.read().splitlines()}

            for compPosition in comparisonData:

                compRef = comparisonData[compPosition][0]  # Should be equal, add sanity check?
                compAlt = comparisonData[compPosition][1]

                if compAlt == 'N':
                    continue

                for pancPosition in pancovData:

                    pancAlt = pancovData[pancPosition][1]

                    pancRef = pancovData[pancPosition][0]


                    if pancAlt == 'N': #N doesn't count
                        continue

                    if compPosition == pancPosition:

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
                                    pileup[int(
                                        pancPosition)] if int(pancPosition) in pileup else 'no pileup available for this position'
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
                            pileup[int(compPosition)] if int(compPosition) in pileup else 'no pileup available for this position'
                        )
                    )
            # Check for new variants (exclusively detected by pancov)
            for pancPosition in pancovData:

                pancAlt = pancovData[pancPosition][1]
                pancRef = pancovData[pancPosition][0]

                if pancAlt == 'N':
                    continue

                for compPosition in comparisonData:
                    if pancPosition == compPosition: #already processed
                        break
                else:
                    #strand bias check 20/80
                    if int(pancPosition) in pileup:
                        if pileup[int(pancPosition)][1] < 0.35 or pileup[int(pancPosition)][1] > 0.65:
                            continue

                    totalNew += 1
                    outfile.write(
                        '{}\t{}\t{}\t{}\n'.format(
                            pancPosition,
                            '{}->{}'.format(pancRef, pancAlt),
                            'Missing',
                            pileup[int(pancPosition)] if int(pancPosition) in pileup else 'no pileup available for this position'
                        )
                    )

    outfile.write(
            'New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew, totalChanged, totalMissed))
