call with --use-singularity --use-conda
report with --report (creates html that collects all results)

sample barcode08 is 'mixed'

sample barcode09 is 'pure'



1. Bauen der Kandidaten-Allel-Liste
2. Für jeden Nanopolish / Medaka-Call in jedem Sample, entscheiden, ob das gecallte Allele HOM ist (kein HET-Call durch Medaka; Coverage; samtools mpileup)
3. Basierend auf allen HOM-Calls in allen Samples die Emission Probabilities für die Allele an ihren jeweiligen Positionen lernen
4. Genotyping über alle Kandidaten-Allele über alle Samples
5. (das klappt nur nicht für Allele, die nur als HET-Calls auftreten, aber nie als HOMs)
6. Ein Sample wird dann als "mixed" gelabelt, wenn es mindestens einen HET-Call gibt
7. Dann noch Ausgabe der Minor Allele Frequencies über alle HET Calls in einem Sample