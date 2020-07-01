import vcf

originalVCF = [r for r in vcf.Reader(open(snakemake.input['originalVCF'], 'r'))]
newVCF = [r for r in vcf.Reader(open(snakemake.input['newVCF'], 'r'))]

def genotypesToText(samples):
    ret = ''
    for sample in samples:
        ret+=sample.gt_bases
    return ret

with open(snakemake.output[0],'w') as outfile:

    outfile.write('{}\t{}\t{}\n'.format('POS','ORIG','NEW'))

    processedPositions = []

    for record in newVCF:
        for originalRecord in originalVCF:
            if record.POS == originalRecord.POS:
                if record.num_het == originalRecord.num_het:
                    #same in both vcfs
                    break
                else:
                    outfile.write('{}\t{}\t{}\n'.format(record.POS,genotypesToText(originalRecord.samples),genotypesToText(record.samples)))
                break
        else:
                outfile.write('{}\t{}\t{}\n'.format(record.POS,'Missing',str(record.alleles[0])+'->'+str(record.alleles[1])))

    for originalRecord in originalVCF:
        for record in newVCF:
            if record.POS == originalRecord.POS:
                break
        else:
                outfile.write('{}\t{}\t{}\n'.format(originalRecord.POS,genotypesToText(originalRecord.samples),'Missing'))
