# Pancov

Pancov is a workflow that builds on the ARTIC-Pipeline (link) to improve nCoV-19 variant calling by aligning reads to a pangenomic variant graph.

# How To

 <a href='https://github.com/esteinig'><img src='doc/logo.png' align="right" height="270" /></a>

This repository is a snakemake workflow that uses conda as a software versioning tool.

In order to run it you need to have a working conda installation as well as snakemake installed.

To execute the workflow, simply edit the ```config.example.yaml``` and rename it to ```config.yaml```.

Then run:

```
snakemake --use-conda
```