
## Install required packages if necessary ----

#packages <- c("dplyr", "here", "readr", "tidyr","diann", "tidyverse")

#biopackgs <- c("MSstats")

#if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
#      install.packages(setdiff(packages, rownames(installed.packages())))  
#}

#if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
#      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#      
#      BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
      
#}


# Modified on 16.02.2021 ----
#install.packages("devtools")
library(tidyverse)
library(dplyr)

##Genes.MaxLFQ gives already normalized MQ-LFQ data for every protein group and is recommended
## Genes.MaxLFQ would be used to extract the quantitative information per feature from the DIANN output
### DIA-NN to MSstats formatting for LFQ normalization ----
## Note: I give full credites to the original author of 
##   this code. Check:https://github.com/MiguelCos/MSstats_labelfree_preprocessing/blob/master/R/diann2msstats.R
#
# I only made a few minor modifications

DIANN_to_MSstats <- function(diann_data, annotation_file){

#remove filepath from File.Name
diann_data1 <- mutate(diann_data, File.Name = str_replace(diann_data[[1]], ".*\\\\", ""))
#diann_data1 <- mutate(diann_data1, File.Name = str_replace(diann_data1[[1]], ".raw.mzml$", ""))
diann_data1 <- mutate(diann_data1, File.Name = str_replace(diann_data1[[1]], ".raw$", ""))

#select relevant columns -> choose which parameter to use for quantitation

for_msstats_prot <- diann_data1 %>% 
      group_by(Stripped.Sequence, 
               Protein.Group, 
               Precursor.Charge, 
               File.Name, 
               Genes.MaxLFQ) %>% 
      summarize()

for_msstats_prot1 <- mutate(for_msstats_prot,
                            File.Name = str_remove(File.Name, ".raw$"))

#merge with MSstats annotation file
annotation <- annotation_file
colnames(annotation)[colnames(annotation) == "Filename"] <- "File.Name"
for_msstats_prot2 <- left_join(for_msstats_prot1, annotation, by = "File.Name")

#change column names to MSstats format
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "BioReplicate"] <- "BioReplicate"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Stripped.Sequence"] <- "PeptideSequence"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Protein.Group"] <- "ProteinName"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Precursor.Charge"] <- "PrecursorCharge"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "File.Name"] <- "Run"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Precursor.Quantity"] <- "Intensity"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Genes.MaxLFQ"] <- "Intensity"

#add missing columns
for_msstats_prot2$IsotopeLabelType <- "L"
for_msstats_prot2$ProductCharge <- NA
for_msstats_prot2$FragmentIon <- NA

#reorder to MSstats format
for_msstats_prot3 <- for_msstats_prot2[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                                           "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                                           "Condition", "BioReplicate", "Run", "Intensity")]

return(for_msstats_prot3)
}


######
#'@ descriptions This is not a true  max quant format converter.
#' To be precise, this is a format similar to max quant output, 
#' Mainly we want to convert diann output to be processed by Perseus
#' So the format is "flexible". We want to make it minimal and extensible
#' We can always added more data to the output.
#'@ param diann_dat a data frame holding the diann output table
#'@ sample_name is the the clean/abbreviated identifier for
#' each sample, which is needed for each row of the data table
#' It could be the same as the run of File.name. Note: this is 
#' NOT the biological sample name, but each physical sample to 
#' run on MS (meaning they could be the technical samples)
DIANN2MQuant <- function(diann_dat, sample_name)
{
  #here, select first sets of fields /columns in the output

  dfm <- diann_dat
  dfm$sample_name = sample_name 
  dfm <- dfm %>% select(sample_name,Run, Protein.Group,  Protein.Ids, Protein.Names, Genes)
)

  return (dfm)
}