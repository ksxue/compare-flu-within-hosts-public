# compare-flu-within-hosts-public
Comparison of within-host diversity across multiple deep-sequencing datasets

This repository contains code and small intermediate data files associated with a comparison of within-host diversity across multiple influenza deep-sequencing datasets. All scripts are meant to be run from the top level of the repository.

Overview
--------

The directory is organized as follows:

    project
    |- README			# Overall description of repository.
	|
	|- analysis			# Custom R, python, and shell scripts to perform analyses in the manuscript.
	|					# Each analysis directory contains a "Run.sh" script that performs all analyses.
    |  |- DownloadDataCallVariants	# Download and organize raw data, map to influenza genomes, and summarize reads at each site.
	|  |- PoonReconstruction		# Map reads from Poon et al. study as single-end reads and analyze diversity.
	|  |- CompareDiversity			# Call variants at a comparable threshold in all datasets and compare their abundance.
	|  |- ReadPairs					# Analyze read pairing in the Poon et al. raw data.
	|  |- CompareR1R2				# Compare genetic diversity on read 1 and read 2.
	|  |- figures					# Generate the figures for the accompanying manuscript.
	|- data				# Small data files, primarily associated with data organization.
	|  |- metadata					# Contains metadata and organizational tracking information for each study.
	|- pipelines		# Standardized computational pipelines for alignment and analysis of sequencing data.
	|- reference		# Reference genomes for H3N2, seasonal H1N1, and pandemic H1N1 influenza.
	|- scripts			# Custom C++, R, and shell scripts for calling and annotating variants.

    
