import os

# Concat

os.system('cat '+ snakemake.input['ref'] + ' ' + snakemake.input['fasta'] + ' > ' + snakemake.output['combined'])

# Muscle
alignmentfile = snakemake.output['realignment']
musclecmd = 'muscle -in ' + snakemake.output['combined'] + ' -clwout ' + alignmentfile
#print(musclecmd)
os.system(musclecmd)

infofile = snakemake.output['info']
with open(alignmentfile, 'r') as infile, open(infofile, 'w') as outfile:
	# Fetch and skip first three lines
	lines = infile.read().splitlines()[3:]
	lastStart = len(lines) - 1 - 2
	currentLine = 0

	position = 0

	while (currentLine <= lastStart):

		reference = lines[currentLine].split()[1]
		variant = lines[currentLine + 1].split()[1]

		for r, v in zip(reference, variant):

			if r != v:
				outfile.write('{}\t{}\t{}'.format(position, r, v))
				if v == 'N':
					outfile.write('\t(Masked)')
				outfile.write('\n')
			if r != '-':
				position += 1

		currentLine += 4


