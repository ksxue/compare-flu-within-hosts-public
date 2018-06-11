require(tidyverse)
require(cowplot)
require(purrr)
require(readr)
require(foreach)
require(seqinr)
require(stringr)

# List output directory.
outdir <- "analysis/CompareDiversity/out/"

# List the column names for variant call files generated directly through
# the variant hard filter script.
# Note that variants that have undergone post-processing will have
# different column names and orders.
VariantCols <- c("Sample","Chr","GenomePos","Gene","Pos","Consensus","Base",
                 "Codon","RefAA","AltAA","Count","AvgQ","AvgReadPos",
                 "Coverage","Freq")

# List the column names for standard samplesheets with reference information.
SamplesheetCols <- c("Raw1","Raw2","Trimmed1","Trimmed2","Sample",
                     "Study","Reference","Subtype")

# Determine the lengths of each gene segment sequenced for each subtype.
GenomeH3N2 <- read.fasta("reference/flu-H3N2/H3N2-Victoria-2011.fasta")
GenomepdmH1N1 <- read.fasta("reference/flu-pdmH1N1/pdmH1N1-California-2009.fasta")



# Dinis et al. 2016 -------------------------------------------------------

# Process variants from the Dinis et al. 2016 study.
# This study sequenced the HA gene from 114 H3N2 and pdmH1N1 patients.

# Based on the 1000 reads mapped in order to infer sample subtype,
# calculate the approximate mapping percentage of the Debbink samples.
DinisMapping <- read.table("analysis/DownloadDataCallVariants/flu-Dinis-subtypes.data",
                             header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(c("Sample","Subtype","MappingRate")) %>%
  mutate(Study="flu-Dinis")

# Determine the total number of H3N2 and pdmH1N1 samples sequenced.
# Note that under the current workflow, only samples with variants
# show in the list of variants.
DinisSamplesheet <- read.table("data/metadata/flu-Dinis/samplesheet-refs.txt",
                               header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% group_by(Subtype) %>%
  summarize(NumSamples=n())

# Read in the Dinis H3N2 and pdmH1N1 variants.
DinisH3N2 <- read.table(paste0(outdir,"flu-Dinis-H3N2.variants"),
                        header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(VariantCols) %>% 
  mutate(Study="flu-Dinis", Subtype="H3N2",
         NumSamplesSequenced=(DinisSamplesheet %>% filter(Subtype=="H3N2"))$NumSamples,
         BpSequenced=length(GenomeH3N2$`4-HA`))

DinispdmH1N1 <- read.table(paste0(outdir,"flu-Dinis-pdmH1N1.variants"),
                        header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(VariantCols) %>% 
  mutate(Study="flu-Dinis", Subtype="pdmH1N1",
         NumSamplesSequenced=(DinisSamplesheet %>% filter(Subtype=="pdmH1N1"))$NumSamples,
         BpSequenced=length(GenomepdmH1N1$`4-HA`))

# Determine the number of variants per sample.
# Make sure to include the samples that have no variants.
# First, obtain the full list of samples sequenced in this study
# and create a dummy dataframe with NumVariants=0 for each sample.
DinisSamples <- read.table("data/metadata/flu-Dinis/samplesheet-refs.txt",
                                     header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% dplyr::select(Sample, Subtype) %>%
  group_by(Subtype) %>% 
  mutate(NumVariants=0, NumSamplesSequenced=n(),
         BpSequenced=ifelse(Subtype=="H3N2",length(GenomeH3N2$`4-HA`),
                            length(GenomepdmH1N1$`4-HA`)))
# Create a list of variants per sample from the list of variants.
# Note that this will not include samples in which no variants were called.
DinisVariantsPerSample <- rbind(
  DinisH3N2 %>% group_by(Sample, Subtype, NumSamplesSequenced, BpSequenced) %>% 
    summarize(NumVariants=n()),
  DinispdmH1N1 %>% group_by(Sample, Subtype, NumSamplesSequenced, BpSequenced) %>%
    summarize(NumVariants=n()))
# Use the anti-join command to determine which samples have no variants.
# Use the union command to join this list of samples to the main list of samples.
DinisVariantsPerSample <- union(DinisVariantsPerSample,
                                anti_join(DinisSamples, DinisVariantsPerSample, 
                                          by=c("Sample"="Sample"))) %>%
  mutate(Study="flu-Dinis")
# Verify that each sample is accounted for uniquely in this count of
# variants per sample.
if((nrow(DinisVariantsPerSample) != nrow(DinisSamples)) |
   (n_distinct(DinisVariantsPerSample$Sample) != 
    n_distinct(DinisVariantsPerSample$Sample))){
  print("Error in tallying variants per sample.")
}

# Poon et al. 2016 --------------------------------------------------------

# Process variants from the Dinis et al. 2016 study.
# This study sequenced the whole flu genome from H3N2 and pdmH1N1 patients.
# These variants have already been pre-processed.
# They consist of averaged frequencies between two sequencing replicates.
# Samples without sequencing replicates were discarded.

# Based on the 1000 reads mapped in order to infer sample subtype,
# calculate the approximate mapping percentage of the Debbink samples.
PoonMapping <- read.table("analysis/DownloadDataCallVariants/flu-Poon-subtypes.data",
                             header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(c("Sample","Subtype","MappingRate")) %>%
  mutate(Study="flu-Poon")

# Determine which samples had sequencing replicates available.
PoonReplicates <- read.table("data/metadata/flu-Poon/Synapse-samplenames-replicates.txt",
                             header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(c("Sample"))

# Determine the total number of H3N2 and pdmH1N1 samples sequenced.
# Note that under the current workflow, only samples with variants
# show in the list of variants.
# Count only samples with two sequencing replicates,
# and count both sequencing replicates as coming from a single sample.
# Disregard one sample for which the two sequencing replicates
# map to different flu subtypes.
PoonSamplesheet <- read.table("data/metadata/flu-Poon/samplesheet-refs.txt",
                               header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% filter(Sample %in% PoonReplicates$Sample) %>%
  separate(Sample, into=c("Household","Visit","Replicate")) %>%
  mutate(Sample=paste(Household,Visit,sep="-")) %>%
  group_by(Sample) %>% filter(n_distinct(Subtype)==1) %>% ungroup() %>%
  group_by(Subtype) %>%
  summarize(NumSamples=n_distinct(Sample))

# Read in the Poon H3N2 and pdmH1N1 variants.
PoonH3N2 <- read.table("analysis/PoonReconstruction/out/ReplicateVariants-H3N2.txt",
                       header=TRUE, stringsAsFactors = FALSE) %>%
  rename(Sample=Strain) %>%
  mutate(Study="flu-Poon", Subtype="H3N2",
         NumSamplesSequenced=(PoonSamplesheet %>% filter(Subtype=="H3N2"))$NumSamples,
         BpSequenced=sum(sapply(GenomeH3N2,length)))

PoonpdmH1N1 <- read.table("analysis/PoonReconstruction/out/ReplicateVariants-pdmH1N1.txt",
                       header=TRUE, stringsAsFactors = FALSE) %>%
  rename(Sample=Strain) %>%
  mutate(Study="flu-Poon", Subtype="pdmH1N1",
         NumSamplesSequenced=(PoonSamplesheet %>% filter(Subtype=="pdmH1N1"))$NumSamples,
         BpSequenced=sum(sapply(GenomepdmH1N1,length)))

# Determine the number of variants per sample.
# Make sure to include the samples that have no variants.
# First, obtain the full list of samples sequenced in duplicate in this study
# and create a dummy dataframe with NumVariants=0 for each sample.
PoonSamples <- read.table("data/metadata/flu-Poon/samplesheet-refs.txt",
                           header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% filter(Sample %in% PoonReplicates$Sample) %>%
  separate(Sample, into=c("Household","Visit","Replicate")) %>%
  mutate(Sample=paste(Household,Visit,sep="-")) %>%
  dplyr::select(Sample, Subtype) %>% distinct() %>%
  group_by(Sample) %>% filter(n_distinct(Subtype)==1) %>% ungroup() %>% 
  group_by(Subtype) %>% 
  mutate(NumVariants=0, NumSamplesSequenced=n(),
         BpSequenced=ifelse(Subtype=="H3N2",sum(sapply(GenomeH3N2,length)),
                            sum(sapply(GenomepdmH1N1,length))))
# Create a list of variants per sample from the list of variants.
# Note that this will not include samples in which no variants were called.
PoonVariantsPerSample <- rbind(
  PoonH3N2 %>% group_by(Sample, Subtype, NumSamplesSequenced, BpSequenced) %>% 
    summarize(NumVariants=n()),
  PoonpdmH1N1 %>% group_by(Sample, Subtype, NumSamplesSequenced, BpSequenced) %>%
    summarize(NumVariants=n()))
# Use the anti-join command to determine which samples have no variants.
# Use the union command to join this list of samples to the main list of samples.
PoonVariantsPerSample <- union(PoonVariantsPerSample,
                                anti_join(PoonSamples, PoonVariantsPerSample, 
                                          by=c("Sample"="Sample"))) %>%
  mutate(Study="flu-Poon")
# Verify that each sample is accounted for uniquely in this count of
# variants per sample.
if((nrow(PoonVariantsPerSample) != nrow(PoonSamples)) |
   (n_distinct(PoonVariantsPerSample$Sample) != 
    n_distinct(PoonVariantsPerSample$Sample))){
  print("Error in tallying variants per sample.")
}



# Debbink et al. 2017 -----------------------------------------------------

# Process variants from the Debbink et al. 2017 study.
# This study sequenced the HA gene from 123 H3N2 patients.

# Based on the 1000 reads mapped in order to infer sample subtype,
# calculate the approximate mapping percentage of the Debbink samples.
DebbinkMapping <- read.table("../within-host-meta/analysis/171129-CallVariants/flu-Debbink-subtypes.data",
                             header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(c("Sample","Subtype","MappingRate")) %>%
  mutate(Study="flu-Debbink")

# Determine the total number of H3N2 and pdmH1N1 samples sequenced.
# Note that under the current workflow, only samples with variants
# show in the list of variants.
DebbinkSamplesheet <- read.table("data/metadata/flu-Debbink/samplesheet-refs.txt",
                               header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% 
  filter(Sample %in% (DebbinkMapping %>% filter(MappingRate>0))$Sample) %>%
  group_by(Subtype) %>%
  summarize(NumSamples=n())

# Read in the Debbink H3N2 variants.
DebbinkH3N2 <- read.table(paste0(outdir,"flu-Debbink-H3N2.variants"),
                        header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(VariantCols) %>% 
  mutate(Study="flu-Debbink", Subtype="H3N2",
         NumSamplesSequenced=(DebbinkSamplesheet %>% filter(Subtype=="H3N2"))$NumSamples,
         BpSequenced=sum(sapply(GenomeH3N2,length)))

# Determine the number of variants per sample.
# Make sure to include the samples that have no variants.
# First, obtain the full list of samples sequenced in this study
# and create a dummy dataframe with NumVariants=0 for each sample.
DebbinkSamples <- read.table("data/metadata/flu-Debbink/samplesheet-refs.txt",
                           header=FALSE, stringsAsFactors = FALSE) %>%
  `colnames<-`(SamplesheetCols) %>% dplyr::select(Sample, Subtype) %>%
  group_by(Subtype) %>% 
  mutate(NumVariants=0, NumSamplesSequenced=n(),
         BpSequenced=ifelse(Subtype=="H3N2",sum(sapply(GenomeH3N2,length)),
                            sum(sapply(GenomepdmH1N1,length))))
#DebbinkSamples$Sample <- as.numeric(DebbinkSamples$Sample)
# Create a list of variants per sample from the list of variants.
# Note that this will not include samples in which no variants were called.
DebbinkVariantsPerSample <- DebbinkH3N2 %>% 
  group_by(Sample, Subtype, NumSamplesSequenced, BpSequenced) %>%
  summarize(NumVariants=n()) %>% ungroup() %>%
  mutate(Sample=ifelse(nchar(Sample)>3,Sample,str_pad(Sample,3,pad="0")))
# Use the anti-join command to determine which samples have no variants.
# Use the union command to join this list of samples to the main list of samples.
DebbinkVariantsPerSample <- union(DebbinkVariantsPerSample,
                                anti_join(DebbinkSamples, DebbinkVariantsPerSample, 
                                          by=c("Sample"="Sample"))) %>%
  mutate(Study="flu-Debbink")
# Verify that each sample is accounted for uniquely in this count of
# variants per sample.
if((nrow(DebbinkVariantsPerSample) != nrow(DebbinkSamples)) |
   (n_distinct(DebbinkVariantsPerSample$Sample) != 
    n_distinct(DebbinkVariantsPerSample$Sample))){
  print("Error in tallying variants per sample.")
}

# Combine all variant files. ----------------------------------------------

Variants <- bind_rows(DinisH3N2, DinispdmH1N1,
                      PoonH3N2, PoonpdmH1N1,
                      DebbinkH3N2)

Mapping <- bind_rows(DinisMapping,
                     PoonMapping,
                     DebbinkMapping)

VariantsPerSample <- bind_rows(DinisVariantsPerSample,
                               PoonVariantsPerSample,
                               DebbinkVariantsPerSample)

# Export lists of variants and mapping rates.
write.table(Variants, paste0(outdir,"variants-all.txt"),
            quote=FALSE, row.names=FALSE)
write.table(Mapping, paste0(outdir,"mapping-all.txt"),
            quote=FALSE, row.names=FALSE)
write.table(VariantsPerSample, paste0(outdir,"variants-per-sample-all.txt"),
            quote=FALSE, row.names=FALSE)

