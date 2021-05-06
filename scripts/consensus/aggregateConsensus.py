from Bio import SeqIO

with open(snakemake.output[0], "w") as outfile:
    records = []

    for string in snakemake.params['input']:
        data = string.split('<::>')
        f = data[0]
        run = data[1]
        barcode = data[2]
        rec = SeqIO.read(f, "fasta")
        rec.id = run + '_' + barcode
        rec.description = ''
        rec.name = ''
        records.append(rec)
    SeqIO.write(records, outfile, "fasta")
