import vcfpy

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

#Step 2 Analysis

def genotypesToText(samples):
    ret = ''
    for sample in samples:
        ret+=str(sample.gt_bases)
    return ret

def altEqual(a1,a2):
    if len(a1) != len(a2):
        return False
    else:
        for x,y in zip(a1,a2):
            if x.serialize() != y.serialize():
                return False
    return True

with open(snakemake.output[0],'w') as outfile:

    totalNew = 0
    totalChanged = 0
    totalMissed = 0

    outfile.write('{}\t{}\t{}\n'.format('Pos', 'Pancov', 'Compared Method'))

    for vcf in mapping:

        outfile.write('// {} / {} //'.format(vcf,mapping[vcf]))

        originalVCF = [r for r in vcfpy.Reader(open(vcf, 'r'))] #pancov
        newVCF = [r for r in vcfpy.Reader(open(mapping[vcf], 'r'))]


        processedPositions = []

        for record in newVCF:
            for originalRecord in originalVCF:
                if record.POS == originalRecord.POS:
                    if altEqual(record.ALT,originalRecord.ALT):
                        #same in both vcfs
                        break
                    else:
                        #changed
                        totalChanged += 1
                        outfile.write('{}\t{}\t{}\n'.format(record.POS,genotypesToText(originalRecord.calls),genotypesToText(record.calls)))
                    break
            else:
                totalMissed += 1
                outfile.write('{}\t{}\t{}\n'.format(record.POS,'Missing',genotypesToText(record.calls)))

        for originalRecord in originalVCF:
            for record in newVCF:
                if record.POS == originalRecord.POS:
                    break
            else:
                totalNew += 1
                outfile.write('{}\t{}\t{}\n'.format(originalRecord.POS,genotypesToText(originalRecord.calls),'Missing'))

    print('New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew,totalChanged,totalMissed))