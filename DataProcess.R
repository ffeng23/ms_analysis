#R code to process the row data
#   -- COPYRIGHT by Feng @ Sensei Bio Jan. 2023
#		--version 1.0
#---
# Tasks: 1) read in the data 
#		 2) preprocess the data
#		 3) reformat to make the format identifical to: i)MaxQuant output data files;
#				ii)Msstate compatible input files
#		 4) write output.
#		 5) (To do in the future) data analysis: clean up, filter, normalize, 
#		 		visualize, statistical comparision.
#---
#


#loading libraries.
#
library(here)
library(dplyr)
library(tidyr)
library(MSstats)
library(MSstatsConvert)
library(diann)

source("functions.R")
#library(readr)

#start reading the files 
df <- diann_load(here("Data","SNS-101_report_short.tsv"))

annotation_file = read.table(file=here("annotation.tsv"),
		header=T, sep="\t")

data.msstats = DIANN_to_MSstats (df, annotation_file)


sample_name <- df$Run %>% 
		sub(pattern="_122[0-9]{1}2022$", replacement="") %>%
		sub(pattern="DIA-",replacement="", fixed=T) %>%
		sub(pattern="CD4\\+T[_]{1,2}iRT_*",replacement="CD4+", fixed=F) %>%
		sub(pattern="00ng",replacement="", fixed=T) %>%
		sub(pattern="SNS-101_VISTA",replacement="SNS_VISTA", fixed=T)  



