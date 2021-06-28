# Pancov

Pancov is a workflow that builds on the ARTIC-Pipeline (link) to improve nCoV-19 variant calling by aligning reads to a pangenomic variant graph.

# Quick Start

<a><img src='doc/logo.png' align="right" height="270" /></a>

This repository is a snakemake workflow that uses conda as a software versioning tool.

In order to run it you need to have a working conda installation as well as snakemake installed.

To execute the workflow, simply edit the ```config.example.yaml``` and rename it to ```config.yaml```.

Place the required input files in ```data/input```.

For example your folder structure could be:

```
data/input/nCoV-2019.reference.fasta
data/input/run1/sample1.medaka.pass.vcf
data/input/run1/sample1.medaka.primertrimmed.rg.sorted.bam
data/input/run1/sample2.medaka.pass.vcf
data/input/run1/sample2.medaka.primertrimmed.rg.sorted.bam
```

Then run:

```
snakemake --use-conda
```

# Methods

For a detailed description of the method refer to: [LINK TO PAPER]
