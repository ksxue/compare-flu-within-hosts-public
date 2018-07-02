#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script organizes raw sequencing files for McCrone et al. 2018.

# Input paths.
BatchSRADownloadFASTQ="pipelines/Batch-SRADownloadFASTQ.sh"
ConsolidateRenameSRA="analysis/DownloadDataCallVariants/ConsolidateRenameSRA-flu-McCrone.sh"
clean=1 # If 0, analysis will run from beginning, overwriting existing intermediates.

# McCrone et al. 2018, eLife
# Download raw FASTQ files from the SRA Bioproject PRJNA412631.
# Use the RunInfo Table for this BioProject to find SRA accessions.
# Rename files to reflect biological sample names.
# Also generate a sample metadata sheet for downstream analyses.

# First, organize the appropriate metadata files.
project="flu-McCrone"
dir="nobackup/${project}/raw"
SRARunTable="data/metadata/${project}/SraRunTable.txt"
numSRA="$(tail -n +2 ${SRARunTable} | wc -l | cut -f1 -d' ')"
runs="data/metadata/${project}/SRR.txt"
mkdir -p nobackup/${project}/raw
mkdir -p nobackup/${project}/trimmed
mkdir -p data/metadata/${project}

# Submit a batch job to download FASTQ files from the SRA.
mkdir -p nobackup/${project}/sge
if [ -z "$(ls ${dir})" ]; then
  echo "Download files from SRA."
  cut -f9 ${SRARunTable} | tail -n +2 > ${runs}
  qsub -cwd -N DownloadMcCrone \
    -t 1-${numSRA} -tc 250 \
    -o nobackup/${project}/sge/ -e nobackup/${project}/sge/ \
    ${BatchSRADownloadFASTQ} ${dir} \
    ${runs} ${clean}
fi

# Rename files to reflect sample names.
# If multiple runs were performed for a single sample,
# then identify those samples and merge the resultant FASTQ files.
if [ "$(ls ${dir}/*.fastq.gz | wc -l)" -eq $(( numSRA * 2 )) ]; then
  
  echo "Rename files and consolidate FASTQs from multiple runs."
  ${ConsolidateRenameSRA} ${SRARunTable} ${dir}

  # Manually combine the remaining exception samples
  # that have three or more corresponding sequencing runs.
  # Then, clear remaining SRR files.
  cat ${dir}/SRR6121221_1.fastq.gz ${dir}/SRR6121245_1.fastq.gz \
    ${dir}/SRR6121244_1.fastq.gz ${dir}/SRR6121243_1.fastq.gz \
	${dir}/SRR6121242_1.fastq.gz > ${dir}/PC1A_1.fastq.gz
  cat ${dir}/SRR6121221_2.fastq.gz ${dir}/SRR6121245_2.fastq.gz \
    ${dir}/SRR6121244_2.fastq.gz ${dir}/SRR6121243_2.fastq.gz \
	${dir}/SRR6121242_2.fastq.gz > ${dir}/PC1A_2.fastq.gz
	
  # Remove all remaining SRR* files.
  rm -f ${dir}/SRR*.fastq.gz
fi

# Download study metadata from associated Github repository.
if [ ! -f data/metadata/${project}/metadata.data ]; then
  wget -O data/metadata/${project}/metadata.data \
    https://raw.githubusercontent.com/lauringlab/Host_level_IAV_evolution/master/data/meta_snv_qual.csv
  # Remove quotes to make file machine-readable.
  sed -i 's/\"//g' data/metadata/${project}/metadata.data
fi

# Prepare dataset sample sheet.
echo "Generate sample data sheet."
if [ ! -f data/metadata/${project}/samplesheet.txt ]; then
  while read sample
  do
    # Pad sample name with zeros if it is a number.
    samplefull=${sample}
    if [[ $sample =~ ^[0-9]+$ ]] ; then
      samplefull=`printf %03d $sample`
    fi
    echo "${dir}/${sample}_1.fastq.gz ${dir}/${sample}_2.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R1.fastq.gz nobackup/${project}/trimmed/${samplefull}_trimmed-R2.fastq.gz ${samplefull} ${project}"
  done < <(cut -f12 ${SRARunTable} | tail -n +2 | sort | uniq ) \
    > data/metadata/${project}/samplesheet.txt
fi



