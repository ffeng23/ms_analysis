
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

  #here, select first sets of fields /columns in the output

  dfm <- diann_dat
  dfm$sample_name = sample_name 
  temp <- dfm %>% group_by(sample_name,Run, Protein.Group,  Protein.Ids, Protein.Names, Genes,
        Stripped.Sequence, First.Protein.Description, Precursor.Charge, Lib.Q.Value, Lib.PG.Q.Value,
        Global.Q.Value, Global.PG.Q.Value,
        
        Genes.MaxLFQ, PG.MaxLFQ, #Q.Value#, #<- this is the quantities will be value to mapped
        PG.Quantity, PG.Normalised, Genes.Quantity, Genes.Normalised
    ) %>% summarise()    ##in this piece of data, we have duplicated with difference in Precursor.Id (different modifications)
                      ## but have same gene quantifications. don't understand!!!
  #switch to wide format  for Genes.MaxLFQ
  temp <- temp %>% pivot_wider(id_cols=c( "Protein.Group",  "Protein.Ids", "Protein.Names", "Genes",
        "Stripped.Sequence", "First.Protein.Description", "Precursor.Charge", "Global.Q.Value", "Global.PG.Q.Value",
        "Lib.Q.Value", "Lib.PG.Q.Value"
        ),
          names_from=sample_name,
          #names_prefix=c("MaxLFQ."),
          values_from=c("Genes.MaxLFQ", "PG.MaxLFQ",#"Q.Value",
            "PG.Quantity", "PG.Normalised", "Genes.Quantity", "Genes.Normalised"),
          )

  #


  return (temp)
}

#'@title remove repeated entries for the same protein/peptide group
#'@descriptions functions to remove duplicated/replicated entries
#'  for each protein and peptide group.
#'@details this function assume the diann format 3 meta columns
#' in the beginning. should be called as the first step after
#'  reading in the data. 
#'  The 'rule' to remove duplicates are defined in the documents 
#'  (at the same folder)
#'@param dat the data frame holding the data of diann format. 
#'  the first 3 columns are "Preotien.Group", Protein.Ids
#'  "Gene"
#''
#'@return a new data frame with no replicated entries
#''
#'
removeReplicates<- function(dat)
{
  #this first step we need to check wether there is gene name for the entry
  # otherwise, we use the protein name
  index.noGeneName<-grep("^\\s*$",x=dat$Genes)
  dat[index.noGeneName,"Genes"]<-dat[index.noGeneName,"Protein.Group"]
  
  #dat is assumed to be of diann format
  # check for the duplicate format
  index.duplicated<-which(duplicated(dat$Genes))

  # get a list of duplicated genes
  genes.duplicated<-unique(dat$Genes[index.duplicated])

  #get a list of 
  list.dup<-list()
  for(i in genes.duplicated)
  {
    #cat("i:",i,"\n")
    list.dup[[i]]<-dat[dat$Genes==i,]
  }

  #now go through the list and do on them to remove duplicates.
  list.dup2<-lapply(list.dup, FUN=removeReplicateEntries)

  dat.new<-dat[!is.element(dat$Genes, genes.duplicated),] 

  dat.new<-rbind(dat.new, do.call(rbind,list.dup2))
  return(dat.new)
}

#this is actual working horse to work on each duplicated group
#remove unnecessary ones and left one entry for each group
# dat1<-dat[dat$Gene=="NACA",]
#dat1<-dat[dat$Gene=="KIF3B",]
# dat1<-dat[dat$Gene=="WASH2P;WASH3P",]
#dat1<-dat[dat$Gene=="POLR1D",]
#MIEF1
#TMPO
#TOR1AIP2
removeReplicateEntries<-function(dat1)
{

  nrow.before=1
  nrow.after=0
  while(nrow.before > nrow.after){
    #first check if all the entries are equal on the measurements
    nrow.before=nrow(dat1)
    if(nrow(unique(dat1[,-c(1:3)]))==1)
    {#just pick the first one and done
      return(dat1[1,])
    }

    #meaning they are not all equal, we pick the one with less missing values
    # counting the NAs
    num.nas<-apply(dat1,1,function(x){sum(is.na(x))}) 

    #pick the one with mininum number of NAs
    min.value<-min(num.nas)

    dat1<-dat1[num.nas==min.value,]
    if(nrow(dat1)==1)
    {
      return(dat1)
    }

    #Try the third rule to see whether there is unique peptide entry
    #looking for things that has one entry (without ; separation)
    index.single<-grep(";",x=dat1$Protein.Ids,fixed=T, invert=TRUE)
    if(length(index.single)>0)
    {
      dat1<-dat1[index.single,]
    } #otherwise don't change it.

    if(nrow(dat1)==1)
    {
      return(dat1)
    }
    nrow.after=nrow(dat1)
  }

  #if we are here, meaning we still can not remove duplicates
  # in this case, we simply pick one and issue an warning.
  warning(paste("Try to remove duplicates but fail for entry gene",
      dat1$Gene,".\n\t We pick a random one as output\n"))
  return(dat1[1,])
}