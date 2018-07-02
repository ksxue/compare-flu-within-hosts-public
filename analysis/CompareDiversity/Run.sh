#!/bin/bash

# Script is meant to be run from the top level of the Github repository.
# Script analyzes within-host variants in the Poon et al. raw data
# with the goal of reproducing Figure 2 from the published analyses.

# Input parameters and directory paths.
analysis="CompareDiversity"
outdir="nobackup/${analysis}/"
mkdir -p ${outdir}
mkdir -p analysis/${analysis}/out

# Output a list of within-host variants for each dataset
# and each flu subtype.
# Where applicable, average variant frequencies calculated
# from independent sequencing replicates.

# Dinis et al. 2016
# This dataset has no sequencing replicates and consists only of variants in HA.
# Therefore, simply concatenate the variant files generated for each sample.
cat nobackup/DownloadDataCallVariants/flu-Dinis/variants/*-H3N2.variants \
  > analysis/${analysis}/out/flu-Dinis-H3N2.variants
cat nobackup/DownloadDataCallVariants/flu-Dinis/variants/*-pdmH1N1.variants \
  > analysis/${analysis}/out/flu-Dinis-pdmH1N1.variants

# Debbink et al. 2017
# This dataset contains sequences from H3N2 samples.
# There are no sequencing replicates available,
# so simply concatenate the variant files generated for each sample.
cat nobackup/DownloadDataCallVariants/flu-Debbink/variants/*-H3N2.variants \
  > analysis/${analysis}/out/flu-Debbink-H3N2.variants

# Poon et al. 2016
# This study includes sequencing replicates for most samples.
# Replicate variants were calculated in the PoonReconstruction analysis folder.

# McCrone et al. 2018
# This study contains sequences from H3N2 and pdmH1N1 samples.
# There are no sequencing replicates available,
# so simply concatenate the variant files generated for each sample.
# Note that some samples are related because they originate from the same patient
# or household.
cat nobackup/DownloadDataCallVariants/flu-McCrone/variants/*-H3N2.variants \
  > analysis/${analysis}/out/flu-McCrone-H3N2.variants
cat nobackup/DownloadDataCallVariants/flu-McCrone/variants/*-pdmH1N1.variants \
  > analysis/${analysis}/out/flu-McCrone-pdmH1N1.variants

# Combine all variant call files, as well as estimates of mapping rates.
Rscript analysis/${analysis}/CombineVariants.R