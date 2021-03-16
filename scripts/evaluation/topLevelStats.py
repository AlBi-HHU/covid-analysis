import pandas

with open(snakemake.input['verification'],'r') as verFile,open(snakemake.input['recovery'],'r') as recFile,open(snakemake.output[0],'w') as outfile:

    vd = verFile.read().splitlines()[-1].split()

    falsePositive = int(vd[1])
    positive = int(vd[3])
    truePositive = positive-falsePositive
    precision =truePositive/positive

    rd = recFile.read().splitlines()[-1].split()

    realPositives = int(rd[3])
    recovered = int(rd[1])
    recall = recovered/realPositives

    f1 = 2*(precision*recall)/(precision+recall)

    het = pandas.read_csv(snakemake.input["het_resume"])
    het = het[het["method"] == "pancov"]
    square_diff = het["square_diff"].mean()

    print(f"{snakemake.config['pangenomeMinCovFactor']},{snakemake.config['pangenomeRVTTPath']},{precision},{recall},{f1},{square_diff}", file=outfile)
