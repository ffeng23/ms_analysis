# R code to do a quick analysis to generate a "raw" data file for the manuscript
#   we do the following on the raw input file from DIANN
#       1. read in using the diann package function
#       2. remove the samples from sns-101 group
#       3. remove the entries for SNS-101 sequences and VISTA protein (gene VSIR)
#       4. get unique sequences or entries for "Protein.Group", "Peptide.Id", "Gene" and "modified sequences"
# --- update 6/28/2023
#***************
#-------------------------------------------
#  below are for the previous comments
#
#===========updated 5/25/2023
# add to remove duplicates and then rerun to generate data 
# results/figures in a new folder 
# 
#-------------------------------
# we do heatmap and pca.
# we assume so far the data have been normalized.
# the data are in ../Data folder

library(here)
library(dplyr)
library(tidyr)
library(diann)
library(readr)
library(FactoMineR)
library(factoextra)
library(ggpubr)

source(here("functions.R"))

# output folder for new results after remove duplicates
# 
output.dir<-"Proj2/removeReps"

#if we want to do without removing duplicates, just comment
# out the line say "df<-removeReplicates(df)"
# and reset the output.dir to below.
#output.dir<-"Proj2"


#read the data
df <- diann_load(here("Proj2","SNS-101_VISTA_CD4+T_iRT_1000ng-DIA_12212022_report.tsv"))

#remove SNS101 treatment group sample entries
df.out<- df %>%
    filter(!grepl(x=Run,pattern="SNS-101")) 
    
#remove the entries with protein start with VISTA and and sns101
df.out<- df.out %>%
    filter(!grepl(x=Protein.Group,pattern="SNS-101-LC"))
df.out<- df.out %>%
    filter(!grepl(x=Genes,pattern="VSIR"))

length(
unique(df.out$Protein.Id))
#[1] 9732
 length(
unique(df.out$Protein.Group))
#[1] 8225

#Genes: 8213
write_tsv(df.out, file=here("Proj2","SNS-101_VISTA_CD4+T_iRT_1000ng-DIA_12212022_report_clean.tsv"))

