#!/bin/bash

# Given the path to a single-end FASTQ file,
# use the FASTQ headers to separate read 1's and
# read 2's into separate files for further analysis.

fastq="$1" # Give path to single-line FASTQ file.
r1="$2" # Give the desired output path for read 1.
r2="$3" # Give the desired output path for read 2.

zcat ${fastq} | paste - - - - -d'\t' | \
  grep /1 | awk '{ OFS="\n"; print $1,$2,$3,$4 }' | gzip -f \
  > ${r1}
zcat ${fastq} | paste - - - - -d'\t' | \
  grep /2 | awk '{ OFS="\n"; print $1,$2,$3,$4 }' | gzip -f \
  > ${r2}
  