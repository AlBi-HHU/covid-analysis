import vcfpy
import os

# Parse other contribs
print('Parsing VCF contributions of all other methods')
others = {}

for othvcf in snakemake.input['other_calls']:
    vcf_reader = vcfpy.Reader(open(othvcf, 'r'))
    for record in vcf_reader:
        if not record.POS in others:
            others[record.POS] = []
        for allele in record.ALT:
            if not allele.value in others[record.POS]:
                others[record.POS].append(allele.value)

finals = {}
# Parsing final variant calls
print('Parsing final VCFs')
for filvcf in snakemake.input['filter_vcfs']:
    vcf_reader = vcfpy.Reader(open(filvcf, 'r'))
    for record in vcf_reader:
        if not record.POS in finals:
            finals[record.POS] = []
        for allele in record.ALT:
            if not allele.value in finals[record.POS]:
                finals[record.POS].append(allele.value)

# Comparing CTX contribs

print('Comparing CTX calls ...')

ctx = {}
unique = {}
used = {}

stats_totals = 0
stats_unique = 0
stats_used = 0
with open(snakemake.output['contrib'],'w') as outfile:
    for ctxvcf in snakemake.input['cortex_calls']:
        if os.stat(ctxvcf).st_size == 0:
            continue
        else:
            print('Non-Empty CTX File: {} \n'.format(ctxvcf))
        vcf_reader = vcfpy.Reader(open(ctxvcf, 'r'))
        for record in vcf_reader:
            if not record.POS in ctx:
                ctx[record.POS] = []
            for allele in record.ALT:
                print(allele)
                if not allele.value in ctx[record.POS]:
                    outfile.write('{}: Detected {} -> {} \n'.format(record.POS,record.REF,allele.value))
                    ctx[record.POS].append(allele.value)
                    stats_total += 1
                    if not record.POS in used:
                        used[record.POS] = []
                    if not allele.value in used[record.POS]:
                        if allele.value in finals[record.POS]:
                            used[record.POS] = True
                            outfile.write('It was actually used! \n')
                            stats_used += 1
                        else:
                            used[record.POS] = False
                            outfile.write('It was never used! \n')
                    if not record.POS in unique:
                        unique[record.POS] = []
                    if not allele.value in unique[record.POS]:
                        if allele.value in others[record.POS]:
                            unique[record.POS] = False
                            outfile.write('It was already known without CTX! \n')
                            stats_unique += 1
                        else:
                            used[record.POS] = True
                            outfile.write('Only CTX found this! \n')
    outfile.write('Total CTX calls: {} \n'.format(stats_totals))
    outfile.write('Unique CTX calls: {} \n'.format(stats_unique))
    outfile.write('Used CTX calls: {} \n'.format(stats_used))
