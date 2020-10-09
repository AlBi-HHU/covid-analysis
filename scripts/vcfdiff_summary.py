import vcf

#Part 1: We want to map the files correctly
mapping = {}
for vcf in snakemake.input['vcfs']:
    mapping[vcf] = 'unknown'
for vcf in snakemake.input['comparisonVcfs']:
    parsedData = vcf.split('/')
    new_run = parsedData[-3]
    new_barcode = parsedData[-2]

    for originalVcf in mapping:
        parsedData = originalVcf.split('/')
        orig_run = parsedData[-2]
        orig_barcode = parsedData[-1].split('.')[0][7:]

        if orig_run == new_run and orig_barcode == new_barcode:
            mapping[vcf] = originalVcf
            print('mapped {} to {}'.format(vcf,originalVcf))
            break
    else:
        print('Could not map anything to: {}'.format(vcf))

for vcf in mapping:
    if mapping[vcf] == 'unknown':
        print('Could not map anything to: {} '.format(vcf))

#Step 2 Analysis

def genotypesToText(samples):
    ret = ''
    for sample in samples:
        ret+=sample.gt_bases
    return ret

with open(snakemake.output[0],'w') as outfile:

	totalNew = 0
	totalChanged = 0
	totalMissed = 0

    for vcf in mapping:

        originalVCF = [r for r in vcf.Reader(open(vcf, 'r'))]
        newVCF = [r for r in vcf.Reader(open(mapping[vcf], 'r'))]

        outfile.write('{}\t{}\t{}\n'.format('POS','ORIG','NEW'))

        processedPositions = []

        for record in newVCF:
            for originalRecord in originalVCF:
                if record.POS == originalRecord.POS:
                    if record.num_het == originalRecord.num_het:
                        #same in both vcfs
                        break
                    else:
	                    #changed
	                    totalChanged += 1
                        outfile.write('{}\t{}\t{}\n'.format(record.POS,genotypesToText(originalRecord.samples),genotypesToText(record.samples)))
                    break
            else:
                totalNew += 1
                outfile.write('{}\t{}\t{}\n'.format(record.POS,'Missing',str(record.alleles[0])+'->'+str(record.alleles[1])))

        for originalRecord in originalVCF:
            for record in newVCF:
                if record.POS == originalRecord.POS:
                    break
            else:
	            totalMissed += 1
                outfile.write('{}\t{}\t{}\n'.format(originalRecord.POS,genotypesToText(originalRecord.samples),'Missing'))

    print('New Variants: {} Changed Variants: {} Missed Variants: {}'.format(totalNew,totalChanged,totalMissed))