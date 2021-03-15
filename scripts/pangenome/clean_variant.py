import vcfpy
import json

##We do two things here: We remove duplicates and also check if the alignments actually support the variants on a nucleotide basis

from collections import defaultdict

def main(in_vcf, supportFile, min_cov, out_vcf):

    nodeSupport = json.load(open(supportFile,'r'))

    reader = vcfpy.Reader.from_path(in_vcf)

    header = reader.header
    
    pos2var = defaultdict(list)
    for record in reader:
        coverage = record.INFO["VCOV"] + record.INFO["RCOV"]
        pos2var[(record.POS, tuple(record.REF), tuple(record.ALT))].append((coverage, record, record.INFO["VCOV"]))

    header.add_filter_line(vcfpy.OrderedDict([('ID', 'LRS'), ('Description', 'Real support below  {} percent'.format(snakemake.config['pagenomeCutoffRealSupport']))]))
    header.add_info_line({"ID": "MULTIPLE", "Type": "Flag", "Number": "1", "Description": "Pangenome found multiple variant at this position"})
    header.add_info_line({"ID": "REALSUPPORT", "Type": "String", "Number": "1", "Description": "Reads that have at least one Match on the Node they support"})
    writer = vcfpy.Writer.from_path(out_vcf, reader.header)

    for key, values in pos2var.items():

        coverage = values[0][0]
        variant = values[0][1]
        vcov = values[0][2]

        if len(values) != 1:
            values.sort(key=lambda x: x[0], reverse=True)
            coverage = values[0][0]
            variant = values[0][1]
            variant.INFO["MULTIPLE"] = True

        #Check for Real Support
        varpath_raw = variant.INFO["VARPATH"]
        nodeIDs = eval(varpath_raw)
        variant.INFO["VARPATH"] = nodeIDs #? hmmm
        supportVals = []
        for nodeID in nodeIDs:
            supportFraction = (nodeSupport[nodeID][0] / vcov)
            supportVals.append(supportFraction)
            if  supportFraction < snakemake.config['pagenomeCutoffRealSupport']:
                variant.add_filter('LRS')

        variant.INFO["REALSUPPORT"] = '/'.join(str(x) for x in supportVals)

        if coverage > min_cov:
            writer.write_record(variant)


if "snakemake" in locals():
    main(snakemake.input['vcf'],snakemake.input['support'],snakemake.config["pangenomeVarMinCov"], snakemake.output[0])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])
