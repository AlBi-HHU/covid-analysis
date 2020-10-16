import os
alignmentfile = snakemake.input['alignment']
infofile = snakemake.output['info']
with open(alignmentfile, 'r') as infile, open(infofile, 'w') as outfile:
	# Fetch and skip first three lines
	lines = infile.read().splitlines()[3:]
	lastStart = len(lines) - 1 - 2
	currentLine = 0

	position = 1

	mode = 'SNP'
	svstring = ''
	svpos = '?'

	while (currentLine <= lastStart):

		reference = lines[currentLine].split()[1]
		variant = lines[currentLine + 1].split()[1]
		
		for r, v in zip(reference, variant):

			if r != v:

				if r == '-' and mode == 'SNP': #insertion starting
					mode = 'INS'
					svstring = v
					svpos = position
				elif r == '-' and mode == 'INS': #insertion enlongation
					mode = 'INS'
					svstring += v
				elif v == '-' and mode == 'SNP': #deletion starting
					mode = 'DEL'
					svstring = r
					svpos = position
				elif v == '-' and mode == 'DEL': #deletion continues
					mode = 'DEL'
					svstring += r
				else:
					if mode == 'DEL': #deletion ended and is followed by SNP
						outfile.write('{}\t{}\t{}'.format(svpos, svstring, '-'))
						outfile.write('\n')
						mode = 'SNP'
					if mode == 'INS': #insertion ended and is followed by SNP
						outfile.write('{}\t{}\t{}'.format(svpos, '-',svstring))
						outfile.write('\n')
						mode = 'SNP'
					outfile.write('{}\t{}\t{}'.format(position, r,v))
					outfile.write('\n')
			else: #reference and variant are identical
				if mode == 'DEL':  # deletion ended
					outfile.write('{}\t{}\t{}'.format(svpos, svstring, '-'))
					outfile.write('\n')
					mode = 'SNP'
				if mode == 'INS':  # insertion ended
					outfile.write('{}\t{}\t{}'.format(svpos, '-', svstring))
					outfile.write('\n')
					mode = 'SNP'

			if r != '-':
				position += 1

		currentLine += 4
