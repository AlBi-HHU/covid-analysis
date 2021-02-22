with open(snakemake.output[0], "w") as outfile:
    total = 0
    recoveredCount = 0
    nanoDropCount = 0
    ambiguous = 0
    for f in snakemake.input:
        outfile.write(f + "\n")
        with open(f, "r") as infile:
            ll = infile.read().splitlines()[1:]
            for l in ll:
                outfile.write(l + "\n")
                d = l.split()

                recovered = d[3]
                altallele = d[2]

                if recovered == "Reject":
                    continue
                elif recovered == "True":
                    recoveredCount += 1
                elif recovered == "Nanopore":
                    nanoDropCount += 1

                if d[2] in ["R", "Y", "S", "W", "K", "M"]:
                    ambiguous += 1
                total += 1

    outfile.write(
        "recovered {} of {} found variants ({} are ambiguous characters, {} nanopore drop of coverage)".format(
            recoveredCount, total, ambiguous, nanoDropCount
        )
    )
