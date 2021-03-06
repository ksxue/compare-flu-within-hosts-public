require(tidyverse)
require(cowplot)
require(purrr)
require(readr)
require(foreach)

# List output directory.
outdir <- "analysis/CompareR1R2/out/"

# Import the list of samples that contain both sequencing replicates.
Replicates <- read.table("data/metadata/flu-Poon/Synapse-samplenames-replicates.txt",
                         header=FALSE, stringsAsFactors = FALSE)
colnames(Replicates) <- c("Sample")

# Read in all *.variant files in the Poon et al. reanalysis variant directory.
# Select only a specific subtype so that the coordinate remain the same.
foreach(subtype=c("H3N2","pdmH1N1"), .combine=rbind) %do% {
  
  # Read in variant files for the appropriate subtype.
  Data <- read.table(paste0("nobackup/CompareR1R2/variants/variants-",
                            subtype,".txt"), header=FALSE,
                     stringsAsFactors = FALSE)
  colnames(Data) <- c("Sample","Chr","GenomePos","Gene","Pos","Consensus","Base",
                      "Codon","RefAA","AltAA","Count","AvgQ","AvgReadPos",
                      "Coverage","Freq")
  
  # Split off the read 1/2 identifier into a separate field.
  Data <- Data %>% 
    separate(Sample, by="-", into=c("Household","Visit","Strain","Read")) %>%
    mutate(Sample=paste(Household,Visit,Strain, sep="-")) %>%
    dplyr::select(-Household, -Visit, -Strain)
  
  
  # Discard variants originating from samples for which
  # only one sequencing replicate was available.
  Data <- Data %>% filter(Sample %in% Replicates$Sample)
  
  # Determine whether variants originated from the same strain,
  # i.e. biological sample. These biological samples contain identifiers like
  # 661-V10 with a three-digit household number and three-digit visit name,
  # according to the original paper.
  Data <- Data %>% separate(Sample, into=c("Household","Visit","Replicate"),
                            remove=FALSE) %>%
    mutate(Strain=paste(Household,Visit,sep="-")) %>%
    dplyr::select(-Household,-Visit)
  
  # Retain only variants that are identified in both sequencing replicates
  # for a particular sample, and on a particular read.
  Data <- Data %>% group_by(Strain, Read, GenomePos) %>%
    filter(n_distinct(Replicate)==2)
  
  # Discard variants called in both sequencing replicates if the
  # consensus base or the variant base differ between samples,
  # i.e. if they do not represent the same variant.
  Data <- Data %>% group_by(Strain, Read, GenomePos) %>%
    filter(n_distinct(Consensus)==1, n_distinct(Base)==1)
  
  # Average the frequencies of the variants calculated in both sequencing replicates
  # to obtain a complete list of variant sites.
  Variants <- Data %>% ungroup() %>%
    group_by(Strain, Read, Chr, GenomePos, Gene, Pos, Consensus,
             Base, Codon, RefAA, AltAA) %>%
    summarize(Freq=mean(Freq))
  
  # Export the replicate variants identified through this analysis.
  write.table(Variants, paste0(outdir,"ReplicateVariants-",subtype,".txt"),
              quote=FALSE, row.names=FALSE, col.names=TRUE)
}

