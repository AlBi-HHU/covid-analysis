import os

#Load Config File
configfile: "config.yaml"


#We detect the runs that are provided as input, alternatively we can use a user-defined subset of runs (this needs to be toggled in the cofig.yaml)
runs = config['runs'] if config['useSubsetOfRuns'] else glob_wildcards('data/input/{run}').run

#For each run we detect the barcodes it contains (individual samples)
barcodes = {}

#It is possible to use a user defined subset here as well
if config['useSubsetOfBarcodes']:
    barcodes = config['barcodes']
    print("Using the following subset of barcodes: {}".format(barcodes))
else:
    for run in runs:
        #We use the medaka output here to determine the existing barcodes, this is arbritrary, any other file could be chosen for this purpose
        barcodes[run] = glob_wildcards('data/input/'+run+'/barcode{barcode}.medaka.'+config['vcf_suffix']).barcode


### REMOVE LATER (THIS SECTION IS ONLY USED FOR EVALUATION OF THE PIPELINE AND IRRELEVANT FOR PRODUCTIVE USE)
gisaidMapping = {}
gisaidMappingInverse = {}

with open('data/input/mappingRunsGisaid.csv', 'r') as infile:
    lines = infile.read().splitlines()
    for l in lines:
        data = l.split()
        run = data[0]
        barcode = data[1]
        if len(data) != 2: #entry in the table
            gisaidID = data[2]
            if not run in gisaidMapping:
                gisaidMapping[run] = {}
            gisaidMapping[run][barcode] = gisaidID
            gisaidMappingInverse[gisaidID] = run+'_'+barcode

def getGisaidFile(run,barcode):
    if run in gisaidMapping:
        if barcode in gisaidMapping[run]:
            return gisaidMapping[run][barcode]
    return None

for run in runs:
    for barcode in barcodes[run]:
        gisaidid = getGisaidFile(run,barcode)

        if gisaidid:
            gisaidpath =  gisaidid  + ".fasta"

            #print(gisaidpath)
            if os.path.isfile(os.path.join('data/input/gisaidseqs',gisaidpath)):
                continue
            else:
                print('No sequence for gisaid id: {} (run: {} bc: {})'.format(gisaidpath,run,barcode))
                del gisaidMapping[run][barcode]
        else:
            print('No mapping for run: {} bc: {}'.format(run,barcode))


illuminaMapping = {}
for run in runs:
    illuminaMapping[run] = set()
    for barcode in barcodes[run]:
        if os.path.isdir(os.path.join('data/input/illumina',(run+'_'+barcode))):
            illuminaMapping[run].add(barcode)

### REMOVE LATER END


# This function assembles all the required output files that serve as a target for the snakemake pipeline
def getInput(wildcards):

    inputList = []

    for run in runs:

        #Mandatory: For each sample we curate a consensus sequence
        inputList += expand('data/output/consensus/'+run+'/{barcode}/consensus.fasta', barcode=barcodes[run])
        #We also create an IGV session file for interactive data exploration
        inputList += expand('data/output/IgvSessions/' + run + '/{barcode}.igv.xml',barcode=barcodes[run])

        #Optional Files: For each sample perform quality control on the reads and output a human-readable .html report
        if config['performQualityControl']:
            inputList += expand('data/output/softClippedSeqs/'+run+'/{barcode}.html',barcode=barcodes[run])

        #Optional Files: For each sample we can annotate our detected variants using SnpEff
        if config['VarAnnotSnpEff']:
            inputList += expand('data/auxiliary/pangenome_vc/'+run+'/{barcode}/filter.annoted.vcf', barcode=barcodes[run])

        ### REMOVE LATER (THIS SECTION IS ONLY USED FOR EVALUATION OF THE PIPELINE AND IRRELEVANT FOR PRODUCTIVE USE)

        inputList += ['data/output/evaluation/comparisonFastaBased/nanopolish.eval']
        inputList += ['data/output/evaluation/comparisonFastaBased/medaka.eval']
        inputList += ['data/output/evaluation/comparisonFastaBased/illumina.eval']
        inputList += ['data/output/evaluation/illumina_verification_pancov.html']
        inputList += ['data/output/evaluation/illumina_verification_medaka.html']
        inputList += ['data/output/evaluation/illumina_verification_nanopolish.html']
        inputList += ['data/output/evaluation/illumina_recovery_pancov.html']
        inputList += ['data/output/evaluation/illumina_recovery_medaka.html']
        inputList += ['data/output/evaluation/illumina_recovery_nanopolish.html']
        #VCF Based
        inputList += ['data/output/evaluation/vcfbased/concordance_pancov.html']
        inputList += ['data/output/evaluation/vcfbased/concordance_medaka.html']
        inputList += ['data/output/evaluation/vcfbased/concordance_nanopolish.html']
        inputList += ['data/output/evaluation/heterozygosity']
        inputList += ['data/output/evaluation/toplevelStats_pancov_'+str(config['pangenomeMinCovFactor'])+'_'+str(config['pangenomeRVTThreshold'])+'.eval']
        for barcode in barcodes[run]:


            if run in illuminaMapping and barcode in illuminaMapping[run]:
                pass
            else:
                #print('skipping run {} barcode {} for illumina comparison as we have no seq yet ...'.format(run,barcode))
                continue

            inputList += [
                'data/auxiliary/illuminaVarCalls/'+run+'_'+barcode+'/ivar.vcf.tsv'
            ]
        ### REMOVE LATER END

    return inputList

#Main rule that aggregates all the targets and is used when they are not specified, see function above for output files that are being created
rule all:
    input:
        getInput

include: 'rules/shared.snk'
include: 'rules/errorCorrection.snk'
include: 'rules/readAnalysis.snk'
include: 'rules/variantAnalysis.snk'
include: 'rules/consensus.snk'
include: 'rules/pangenome.snk'
include: 'rules/discovery.snk'
include: 'rules/pangenome_variant_call.snk'
### REMOVE LATER (THE INCLUDES BELOW ARE ONLY FOR EVALUATION)
include: 'rules/pangenome_eval.snk'
include: 'rules/illumina.snk'
