# How To



get data from luke (now too big to include in git)

move data into input folder in following structure

data/input/{run}/result_hac



(rename barcodeX.medaka.primertrimmed.medaka.Y into barcodeX.medaka.primertrimmed.Y to have coherence with nanopolish and avoid "double" indexing)



call with

​    snakemake --use-conda

report with 

​    snakemake --report (creates html that collects all results)





# Various Notes



1. Bauen der Kandidaten-Allel-Liste
2. Für jeden Nanopolish / Medaka-Call in jedem Sample, entscheiden, ob das gecallte Allele HOM ist (kein HET-Call durch Medaka; Coverage; samtools mpileup)
3. Basierend auf allen HOM-Calls in allen Samples die Emission Probabilities für die Allele an ihren jeweiligen Positionen lernen
4. Genotyping über alle Kandidaten-Allele über alle Samples
5. (das klappt nur nicht für Allele, die nur als HET-Calls auftreten, aber nie als HOMs)
6. Ein Sample wird dann als "mixed" gelabelt, wenn es mindestens einen HET-Call gibt
7. Dann noch Ausgabe der Minor Allele Frequencies über alle HET Calls in einem Sample





## interesting positions: (not sure what to do here)

11182 T,C (bc9/24bc)

14408 T,C (bc9/24bc)

​                  (bc15/24bc)

11083 (T variant with DELs)

13225 (no strand bias, non-called variant? bc20/24bc)

20887 (bc 17/24bc)

23683 (bc 18/24bc)





```
outputList:
run/barcode hom/het addedToDBID    frequencies if hom
```