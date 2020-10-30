from shared import *
import statistics

iterator = iter(snakemake.input['iteratorList'])


with open(snakemake.output[0],'w') as outfile:

    totalChanged = 0
    totalMissed = 0
    totalNew = 0
    totalDetectedA = 0
    totalDetectedB = 0
    identicalVars = 0

    meanPileups = {}

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
                totalDetectedB += 1

                for pancPosition in pancovData:


                    pancAlt = pancovData[pancPosition][1]

                    pancRef = pancovData[pancPosition][0]


                    if pancAlt == 'N': #N doesn't count
                        continue

                    if compPosition == pancPosition:

                        if pancAlt == compAlt:
                            identicalVars += 1
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
                totalDetectedA += 1

                for compPosition in comparisonData:
                    if pancPosition == compPosition: #already processed

                        compRef = comparisonData[compPosition][0]  # Should be equal, add sanity check?
                        compAlt = comparisonData[compPosition][1]

                        if compAlt == 'N':
                            continue
                            
                        break
                else:
                    '''
                    #strand bias check 20/80
                    if int(pancPosition) in pileup:
                        strandbias = pileup[int(pancPosition)][1]
                        breakout = False
                        for bias in strandbias:
                            if bias > 80:
                                breakout = True
                                break
                        if breakout:
                            break
                    '''

                    totalNew += 1

                    #check cohort frequencies on demand
                    if not pancPosition in meanPileups:
                        meanPileups[pancPosition] = {}
                    if not pancAlt in meanPileups[pancPosition]:
                        print('fetching median allele frequencies for {} : {}'.format(pancPosition,pancAlt))
                        median = {
                            pancAlt.lower() : [],
                            pancAlt.upper() : []
                        }

                        #Samples that cover the position
                        coveringSamples = 0

                        for f in snakemake.input['pileupOnly']:

                            pileup = parsePileupStrandAware(f)[0]

                            if pancPosition in pileup:
                                coveringSamples += 1
                                if pancAlt.lower() in pileup:
                                    median[pancAlt.lower()] += pileup[pancAlt.lower()]
                                if pancAlt.upper() in pileup:
                                    median[pancAlt.upper()] += pileup[pancAlt.upper()]

                        median = {
                            k : statistics.median(v) for k,v in median.items()
                        }

                        meanPileups[pancPosition][pancAlt] = median

                    outfile.write(
                        '{}\t{}\t{}\t{}\n'.format(
                            pancPosition,
                            '{} -> {} [{}]'.format(pancRef, pancAlt, meanPileups[pancPosition][pancAlt]),
                            'Missing',
                            pileup[int(pancPosition)] if int(pancPosition) in pileup else 'no pileup available for this position'
                        )
                    )
    outfile.write(
            'Total Vars > Pancov: {} Compared Method: {}\n'.format(totalDetectedA,totalDetectedB))
    outfile.write(
            'New Variants: {} Changed Variants: {} Missed Variants: {} Identical Variants: {}'.format(totalNew, totalChanged, totalMissed, identicalVars))
