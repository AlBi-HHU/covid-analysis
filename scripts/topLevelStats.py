with open(snakemake.input['verification'],'r') as verFile,open(snakemake.input['recovery'],'r') as recFile,open(snakemake.output[0],'w') as ouftile:

	vd = verFile.read().splitlines()[-1].split()

	falsePositive = int(vd[1])
	positive = int(vd[3])
	truePositive = positive-falsePositive
	precision =truePositive/positive

	rd = recFile.read().splitlines()[-1].split()

	realPositives = int(rd[1])
	recovered = int(rd[3])
	recall = recovered/realPositives

	f1 = 2*(precision*recall)/(precision+recall)

	snakemake.output('step:{},precision:{},recall:{},f1:{}'.format(snakemake.config['pangenomeMinCovFactor'],precision,recall,f1))

