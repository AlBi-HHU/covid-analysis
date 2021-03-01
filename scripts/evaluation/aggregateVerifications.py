
with open(snakemake.output[0],'w') as outfile:
    total = 0
    rejected = 0
    filteredIllumina = 0
    filteredNanopore = 0
    for f in snakemake.input:
        outfile.write(f+'\n')
        with open(f,'r') as infile:
            ll = infile.read().splitlines()[1:]
            for l in ll:
                outfile.write(l+'\n')
                d = l.split()
                reject = d[3]
                if reject == "True":
                    rejected += 1
                elif reject == "False":
                    total += 1
                elif reject == "Illumina":
                    filteredIllumina += 1
                elif reject == "Nanopore":
                    filteredNanopore += 1

    outfile.write('rejected {} of {} found variants, no decision on {} ({} illumina {} nanopore) variants due to low coverage'.format(rejected, total, filteredIllumina + filteredNanopore, filteredIllumina, filteredNanopore))
