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

#library(readr)

#start reading the files 
df <- diann_load(here("Data","SNS-101_report_short.tsv"))



