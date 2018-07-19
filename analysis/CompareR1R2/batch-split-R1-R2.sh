#!/bin/bash

# This script is meant to be run from the top level of the Github repository.
# Script for batch submission of jobs to split read pairs.

# Input parameters.
pipeline="analysis/CompareR1R2/split-R1-R2.sh"
samplesheet="$1" # Give the path to a samplesheet listing rawread1, rawread2, trimmed1, trimmed2, samplename.
dir="$2" # Give the desired working directory for the read 1 and read 2 files.

# Extract the appropriate line of the list of run IDs to submit to the script.
fastq1=`cut -f1 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract path of the FASTQ file.
sample=`cut -f5 -d' ' ${samplesheet} | awk "NR==$SGE_TASK_ID"` # Extract sample name.
echo ${fastq1}
echo ${sample}
${pipeline} ${fastq1} ${dir}/${sample}-R1.fastq.gz ${dir}/${sample}-R2.fastq.gz