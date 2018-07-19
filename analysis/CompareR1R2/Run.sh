#!/bin/bash

# This script compares the results of analyzing only the first or second reads
# of samples from the Poon et al. dataset.

# Input parameters and directory paths.
project="flu-Poon"
analysis="CompareR1R2"
outdir="nobackup/${analysis}"
samplesheet="data/metadata/flu-Poon/samplesheet-SingleEnd.txt"
numsamples="$(wc -l ${samplesheet} | cut -f1 -d' ')"
BatchInferSubtypeSingleEnd="scripts/InferSubtype/Batch-InferSubtype-SingleEnd.sh"
BatchAlignSummarizeRealignSingleEnd="pipelines/Batch-AlignSummarizeRealign-SingleEnd.sh"
BatchCallVariantsHardFilter="scripts/CallVariants/Batch-CallVariants-HardFilter.sh"
mkdir -p ${outdir}
mkdir -p ${outdir}/raw

# Split FASTQ files into files containing only read 1
# or read 2 reads
if [ ! "$(ls ${outdir}/raw/*.fastq.gz | wc -l)" -eq $(( 2 * ${numsamples} )) ]; then
  echo "Split read 1 and read 2."
  # Submit batch jobs to split read 1 and read 2 for each sample.
  qsub -cwd -N SplitReads \
    -t 1-${numsamples} -tc 250 \
    -o nobackup/sge -e nobackup/sge \
    analysis/${analysis}/batch-split-R1-R2.sh \
    ${samplesheet} ${outdir}/raw
fi

# Generate new samplesheets for the read1 and read2 files
echo "Generate samplesheets."
samplesheetR1R2="data/metadata/flu-Poon/samplesheet-R1R2.txt"
while read fastq1 fastq2 trimmed1 trimmed2 sample study
do
  echo ${outdir}/raw/${sample}-R1.fastq.gz ${outdir}/raw/${sample}-R1.fastq.gz ${outdir}/raw/${sample}-R1.fastq.gz ${outdir}/raw/${sample}-R1.fastq.gz ${sample}-R1 ${project}
  echo ${outdir}/raw/${sample}-R2.fastq.gz ${outdir}/raw/${sample}-R2.fastq.gz ${outdir}/raw/${sample}-R2.fastq.gz ${outdir}/raw/${sample}-R2.fastq.gz ${sample}-R2 ${project}
done < ${samplesheet} > ${samplesheetR1R2}

# Infer the subtype of each flu sample.
dir=${outdir}
if [ ! "$(ls ${dir}/*.map | wc -l)" -eq $(( 2 * ${numsamples} )) ]; then
  echo "Inferring sample subtypes."
  # Submit batch jobs to analyze the subtypes of each flu sample.
  qsub -cwd -N InferSubtype -l h_rt=1:00:00 \
    -hold_jid InferSubtype \
    -t 1-$(( 2 * ${numsamples} )) -tc 500 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchInferSubtypeSingleEnd} ${samplesheetR1R2} ${dir} \
    scripts/InferSubtype/flu-subtypes.data 1000
fi

# Generate a samplesheet with reference genomes.
if [ "$(ls ${dir}/*.map | wc -l)" -eq $(( 2 * ${numsamples} )) ] && \
  [ ! -f data/metadata/${project}/samplesheet-R1R2-refs.txt ]; then
  echo "Generate sample sheet with inferred subtypes."
  # Generate a summary of the best-matching reference genomes
  # for each sample and the mapping rates to these genomes.
  while read fastq1 fastq2 trimmed1 trimmed2 sample study
  do
    head -n 1 ${dir}/${sample}.map
  done < ${samplesheetR1R2} > ${dir}/${project}-R1R2-subtypes.summary
  cut -f3,5,6 -d' ' ${dir}/${project}-R1R2-subtypes.summary \
    > analysis/${analysis}/${project}-R1R2-subtypes.data

  # Merge information from the samplesheet and about sample subtypes.
  while read fastq1 fastq2 trimmed1 trimmed2 sample study \
    fastq12 fastq22 sample2 reference refname maprate
  do
    if [ "$sample" == "$sample2" ]
    then
      echo ${fastq1} ${fastq2} ${trimmed1} ${trimmed2} \
	    ${sample} ${study} ${reference} ${refname}
    else
      echo "Error in matching subtype to sample."
    fi
  done < <( paste -d' ' ${samplesheetR1R2} ${dir}/${project}-R1R2-subtypes.summary ) \
    > ${dir}/${project}-R1R2.samplesheet
  
	
  # Copy samplesheet with reference genomes to metadata folder.
  cp ${dir}/${project}-R1R2.samplesheet \
    data/metadata/${project}/samplesheet-R1R2-refs.txt
fi

# Align reads to reference genome, infer sample consensus,
# re-align reads to consensus genome, and tally and annotate variants.
if [ ! "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq $(( 2 * ${numsamples} )) ] && \
  [ -f data/metadata/${project}/samplesheet-R1R2-refs.txt ]; then
  echo "Aligning reads and calling variants."
  # Submit batch jobs to map all reads for each flu samples.
  # Map reads and call variants based on initial mapping to the
  # best-match reference genome determined above.
  qsub -cwd -N AlignReads -l h_rt=48:00:00 \
    -t 1-$(( 2 * ${numsamples} )) -tc 500 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchAlignSummarizeRealignSingleEnd} \
    ${dir}/${project}-R1R2.samplesheet ${dir} 1
fi

# Call variants using a hard filter.
# Classify a site as variable if there are at least 200 reads at that site
# and at least 3% of them support a non-consensus base
if [ "$(ls ${dir}/variants/*.variants | wc -l)" -eq 0 ] && \
  [ "$(ls ${dir}/*-realigned-annotated.summary | wc -l)" -eq $(( 2 * ${numsamples} )) ] && \
  [ -f data/metadata/${project}/samplesheet-R1R2-refs.txt ]; then
  echo "Calling variants."
  mkdir -p ${dir}/variants
  # Submit batch jobs to call variants in each sample individually.
  qsub -cwd -N CallVariantsHardFilter -l h_rt=48:00:00 \
    -t 1-$(( 2 * ${numsamples} )) -tc 500 \
    -o nobackup/sge -e nobackup/sge \
    ${BatchCallVariantsHardFilter} \
    data/metadata/${project}/samplesheet-R1R2-refs.txt ${dir} \
    ${dir}/variants 0.03 200 
fi

# Determine which samples do not produce variant files,
# i.e. presumably have no variable sites under the given criteria.
if [ ! "$(ls ${dir}/variants/*.variants | wc -l)" -eq $(( 2 * ${numsamples} )) ] && \
  [ "$(ls ${dir}/variants/*.variants | wc -l)" > 0 ] && \
  [ ! -f ${dir}/variants/novariants.txt ] && \
  [ -f data/metadata/${project}/samplesheet-R1R2-refs.txt ]; then
  echo "Identifying samples with no variable sites."
  while read fastq1 fastq1 trimmed1 trimmed2 sample project refpath ref
  do
    if [ ! -f ${dir}/variants/${sample}-${ref}.variants ]; then
	  echo ${sample}
	fi
  done < data/metadata/${project}/samplesheet-R1R2-refs.txt \
    > ${dir}/variants/novariants.txt
fi

# Extract variant frequencies in order to reconstruct Figure 2.
Rscript analysis/CompareR1R2/ReconstructFigure2-R1R2.R

# Summarize the variant frequency calls.
cat nobackup/${analysis}/variants/*H3N2* \
  > nobackup/${analysis}/variants/variants-H3N2.txt
cat nobackup/${analysis}/variants/*pdmH1N1* \
  > nobackup/${analysis}/variants/variants-pdmH1N1.txt
  
# Extract replicate variants for each sample
# for read 1 alone and read 2 alone.
Rscript analysis/${analysis}/CallReplicateVariants-flu-Poon-R1R2.R
Rscript analysis/${analysis}/CombineVariants-R1R2.R