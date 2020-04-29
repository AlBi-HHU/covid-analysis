#9 -> clean 8 -> mixed

configfile: "config.yaml"

ks = config['ks']

methods = ['medaka','nanopolish']
#Uses the IGV sessions which is completely arbitrary, could use any other input file here to get the barcode ids

runs = config['runs'] if config['useSubsetOfRuns'] else glob_wildcards('data/input/{run}/result_hac').run

barcodes = {}

if config['useSubsetOfBarcodes']:
    barcodes = config['barcodes']
else:
    for run in runs:
        barcodes[run] = glob_wildcards('data/input/'+run+'/result_hac/barcode{barcode}.medaka.primertrimmed.vcf').barcode


def getInput():
    inputList = []

    inputList.append([
                      'data/output/tobigram.svg',
                      'data/auxiliary/interestingHOMPositions.json',
                      'data/auxiliary/interestingPositions.json',
                      'data/auxiliary/sampleClassification.tsv'
                      ])

    for run in runs:
        if config['generateGFAs']:
            inputList += expand('data/auxiliary/graphs/{method}/'+run+'/{barcode}/{k}.gfa',method=methods,barcode=barcodes[run],k=ks)

        if config['performQualityControl']:
            inputList += expand('data/output/softClippedSeqs/{method}/'+run+'/{barcode}.html',method=methods,barcode=barcodes[run])

        if config['generateKmerProfiles']:
            inputList += expand('data/output/kmerHistograms/{method}/'+run+'/{barcode}_{k}.svg',method=methods,barcode=barcodes[run],k=ks)

        if config['performCorrections']:
            inputList += [
                          expand('data/output/{corrections_type}/{method}/'+run+'/{k}/{barcode}.svg',
                                 method=methods,
                                 barcode=barcodes[run],
                                 k=ks,
                                 corrections_type=["corrections", "corrections_clip"])
                          for ks in range(
                              config["correctionMinKmer"],
                              config["correctionMaxKmer"] + 1,
                              config["correctionKmerStep"]
                              )
                            ]

        if config['generatePangenome']:
           inputList += ["data/auxiliary/pangenome/all_reads_and_ref.paths.gfa"]

        if config['pangenomeVariantCalling']:
            inputList += expand('data/auxiliary/pangenome_vc/{method}/'+run+'/{barcode}/variant.vcf', method=methods, barcode=barcodes[run])

        if config['realignment']:
            inputList += expand('data/auxiliary/realignment/{corrections_type}/{method}/'+run+'/{barcode}.sorted.bam.bai', corrections_type=["corrections", "corrections_clip"], method=methods, barcode=barcodes[run])

        inputList += expand('data/output/IgvSessions/{method}/'+run+'/{barcode}.igv.xml',method=methods,barcode=barcodes[run])

    return inputList


rule all:
    input:
        getInput()


include: 'rules/errorCorrection.snk'
include: 'rules/debruijn.snk'
include: 'rules/kmerAnalysis.snk'
include: 'rules/variantAnalysis.snk'
include: 'rules/pangenome.snk'
include: 'rules/realignment.snk'
include: 'rules/pangenome_variant_call.snk'
