import os


with open('comparisonFiles.txt','r') as comparisonFilesFile:
	lines = comparisonFilesFile.read().splitlines()
	for line in lines:
		data = line.split('\t')
		run = data[0]
		barcode = data[1]

		os.system(
			'snakemake --use-conda --cores all --config config.yaml \
			--config useSubsetOfRuns=True \
			--config runs=[{}] \
			--config useSubsetOfBarcodes=True \
			--config barcodes= \{ {} : [{}] \} \
			-R concatenate_vcf \
			data/auxiliary/pangenome_vc/{}/{}/filter.vcf'.format(
				run,
				run,
				barcode,
				run,
				barcode
			)
		)