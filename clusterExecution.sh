#!/bin/bash

module load Miniconda/3_snakemake
module load Snakemake/5.10.0

while getopts "p:j:" opt #r:
do
   case $opt in
       p) projectID="$OPTARG" ;;
       j) maxNrOfConcurrentJobs="$OPTARG" ;;
       #r) reportFile="$OPTARG" ;;
   esac
done

if [ -z "$projectID" ] || [ -z "$maxNrOfConcurrentJobs" ] #|| [ -z "$reportFile" ]
then
   echo "Usage: clusterExecution -p projectID -j maxNrOfConcurrentJobs -r reportFile"
   exit 1
fi

mkdir -p clusterLogs

type snakemake >/dev/null 2>&1 || { echo >&2 "I require snakemake but it's not installed or added to your path.  Aborting..."; exit 1; }

snakemake --jobs $maxNrOfConcurrentJobs --use-conda --reason --jobscript cluster/jobscript.sh --cluster "qsub -e clusterLogs/{rule}.{wildcards}.{jobid}.errors -o clusterLogs/{rule}.{wildcards}.{jobid}.output -A ${projectID} -l select=1:ncpus={params.cpus}:ngpus={params.gpus}:mem={params.mem} -l walltime={params.walltime}"

#snakemake --report $reportFile
