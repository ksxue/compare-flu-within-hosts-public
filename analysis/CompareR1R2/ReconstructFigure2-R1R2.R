require(tidyverse)
require(cowplot)
require(purrr)
require(readr)
require(foreach)

# Reconstruct Figure 2 from the published Poon et al. paper as closely as possible.
# In this analysis, use only the read 1's or only the read 2's to reconstruct the figure.
# To do this, first identify variable sites in the H3N2 samples that are identified
# at a frequency above 0.03 at sites with coverage above 200 in both replicates
# in the full analysis (i.e. with both read 1 and read 2 reads).
# Extract the full read count information for all H3N2 samples at those variable sites
# for read 1 alone and read 2 alone.
# Retain only information for samples with both sequencing replicates available.

# List output directory.
outdir <- "analysis/CompareR1R2/out/"

# Read in the list of variants identified in both replicates of H3N2 influenza samples.
Variants <- read.table("analysis/PoonReconstruction/out/ReplicateVariants-H3N2.txt",
                       header=TRUE, stringsAsFactors = FALSE)

# Retain only HA variants that are observed in more than one
# distinct biological sample.
Variants <- Variants %>% group_by(GenomePos) %>%
  filter(n_distinct(Strain)>1) %>%
  filter(Gene=="4-HA", Codon)

# Read in all of the raw read count files for the H3N2 variants.
# Do this for both read 1 alone and read 2 alone.
subtype <- "H3N2"
path <- "nobackup/CompareR1R2/"
files <- dir(path, pattern=paste0(subtype,"-realigned-annotated.summary$"))

# Combine all of the variant files into a single dataframe.
Data <- files %>%
  map(~ read.table(file.path(path, .), 
                   header=FALSE, stringsAsFactors = FALSE)) %>%
  reduce(rbind)
colnames(Data) <- c("Chr","Pos","Base","RefBase","GenomePos","Count","AvgQ","AvgReadPos",
                    "Gene","Codon","RefAA","AltAA","Syn","FourfoldSyn",
                    "Sample")
# Split off the read 1/2 identifier into a separate field.
Data <- Data %>% 
  separate(Sample, by="-", into=c("Household","Visit","Strain","Read")) %>%
  mutate(Sample=paste(Household,Visit,Strain, sep="-")) %>%
  dplyr::select(-Household, -Visit, -Strain)
Data <- Data %>% filter(Gene != "none") %>% group_by(Sample, Read, Gene, Pos) %>%
  mutate(Consensus=Base[which.max(Count)], Freq=Count/sum(Count),
         Coverage=sum(Count))

# Retain only sites of HA variation in two or more samples,
# as determined previously.
DataVariable <- Data %>% filter(GenomePos %in% unique(Variants$GenomePos))

# Retain only samples for which both sequencing replicates are available.
# Import the list of samples that contain both sequencing replicates.
Replicates <- read.table("data/metadata/flu-Poon/Synapse-samplenames-replicates.txt",
                         header=FALSE, stringsAsFactors = FALSE)
colnames(Replicates) <- c("Sample")
DataVariable <- DataVariable %>% filter(Sample %in% Replicates$Sample)

# For each site and each strain, average the variant frequencies
# across the two sequencing replicates.
DataVariable <- DataVariable %>% 
  separate(Sample, into=c("Household","Visit","Replicate")) %>%
  mutate(Strain=paste(Household,Visit,sep="-")) %>%
  dplyr::select(-Household,-Visit)
DataVariable <- DataVariable %>% ungroup() %>%
  group_by(Strain, Read, Chr, Pos, Base, GenomePos,
           Gene, Codon) %>% 
  summarize(Freq=mean(Freq))

# Export the averaged base frequencies at each variable HA site.
write.table(DataVariable, paste0(outdir,"Figure2Frequencies.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE)
