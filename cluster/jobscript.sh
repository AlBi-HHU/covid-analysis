#!/bin/bash

#Load Snakemake and Singularity for cluster execution if a module system exists on the HPC used for computation
module load Miniconda/3_snakemake
module load Snakemake/5.10.0

#properties = {properties}

{exec_job}
