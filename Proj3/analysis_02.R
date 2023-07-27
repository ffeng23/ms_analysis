# R code to do analysis of MS data with T cell in vitro
#
#----- 7/26/23 ---
#  copied from analysis_01.R to do the second run of analysis
#	VISTA modification data.
#     the data set are two software (DIANN vs. Spectronaut)
#		two condifitions H only amino acid mod vs 6 aa mods.
#	data are at ./Data folder
#
#
#####################################
#   below are the comments from analysis_01.R
#----- Started 7/12/2023
#- in this module we will run analysis on VISTA modifications
#	so that compare ratio of modification at each position for each peptide
##############################

library(here)
library(dplyr)
library(tidyr)
library(diann)
library(readr)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(stringr)

source(here("functions.R"))
source(here("Proj3", "functions_modifications.R"))
# output folder for new results after remove duplicates
# 
output.dir<-"Proj3"
input.dir<-"Proj3/Data"
out.file0<-"VISTA_modifications_stats"

#read the data
df.diann.H <- diann_load(here(input.dir,"DIA-NN_VISTA_Modifications_DEPC(H)_report.tsv"))

df2<- df.diann.H %>% dplyr::select("Run", "Protein.Group","Protein.Ids","Protein.Names","Genes", "PG.MaxLFQ",
		"Genes.MaxLFQ", "Modified.Sequence", "Stripped.Sequence", "Precursor.Id","Precursor.Normalised", "Ms1.Area" )

#parse out the condition/grouping info from the field of Run
df2$Run2<-sub(df2$Run, pattern="_[0-9]+$",replacement="") #<-get rid of date, trailing 
df2$Repeat<-substr(df2$Run2, nchar(df2$Run2), nchar(df2$Run2)) #<- get repeat number
df2$Run2<-sub(df2$Run2,pattern="\\-[1-3]{1}$", replacement="") #<- get rid of repeat number

df2$Run2<-sub(df2$Run2,pattern="^DIA_iRT_", replacement="") #<- get rid of leading "DIA_iRT"



#get rid of Unimod:4 which we don't want to take care of in this round
df2$Modified.Sequence <- gsub(pattern="\\(UniMod\\:[0-9]+\\)",
			x=df2$Modified.Sequence,"")

##
# first we need to go through the modified sequences to get
# all possible modifications
mseq<-str_extract_all(pattern="\\([a-zA-Z0-9]+\\)",
			string=df2$Modified.Sequence )

mods<-unique(unlist(mseq))

#we know now it is two. we get rid of the first one and do the second one and 
#m_remove<-mods[2]
m<-mods[1]
out.file<-paste0(out.file0, "_",m,"_run2.csv")
df2$Modified.Sequence2<-df2$Modified.Sequence
		#get rid of mod[1] "(DEPC)", since we want to study "Fumi"
		# don't need to do this now.
		#  Don't RUN
			df2$Modified.Sequence2<-gsub(pattern=m_remove,
						x=df2$Modified.Sequence,"", fixed=T)
		##################end of don't run section###########


#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]
runModStats(dt=df2, output.dir, out.file, m=m, quantity.field="Ms1.Area")


########################################################3
#now let's the six aa modifications for DEPC
#########################################################

#read the data

df.diann.AA6 <- diann_load(here(input.dir,"DIA-NN_VISTA_Modifications_DEPC(HCKSTY)_report.tsv"))

df2<- df.diann.AA6 %>% dplyr::select("Run", "Protein.Group","Protein.Ids","Protein.Names","Genes", "PG.MaxLFQ",
		"Genes.MaxLFQ", "Modified.Sequence", "Stripped.Sequence", "Precursor.Id","Precursor.Normalised", "Ms1.Area" )

#parse out the condition/grouping info from the field of Run
df2$Run2<-sub(df2$Run, pattern="_[0-9]+$",replacement="") #<-get rid of date, trailing 
df2$Repeat<-substr(df2$Run2, nchar(df2$Run2), nchar(df2$Run2)) #<- get repeat number
df2$Run2<-sub(df2$Run2,pattern="\\-[1-3]{1}$", replacement="") #<- get rid of repeat number

df2$Run2<-sub(df2$Run2,pattern="^DIA_iRT_", replacement="") #<- get rid of leading "DIA_iRT"



#get rid of Unimod:4 which we don't want to take care of in this round
df2$Modified.Sequence <- gsub(pattern="\\(UniMod\\:[0-9]+\\)",
			x=df2$Modified.Sequence,"")

##
# first we need to go through the modified sequences to get
# all possible modifications
mseq<-str_extract_all(pattern="\\([a-zA-Z0-9]+\\)",
			string=df2$Modified.Sequence )

mods<-unique(unlist(mseq))

#we know now it is two. we get rid of the first one and do the second one and 
#m_remove<-mods[2]
m<-mods[1]
out.file<-paste0(out.file0, "_",m,"_aa6_run2.csv")
df2$Modified.Sequence2<-df2$Modified.Sequence
		#get rid of mod[1] "(DEPC)", since we want to study "Fumi"
		# don't need to do this now.
		#  Don't RUN
			df2$Modified.Sequence2<-gsub(pattern=m_remove,
						x=df2$Modified.Sequence,"", fixed=T)
		##################end of don't run section###########
#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]
runModStats(dt=df2, output.dir, out.file,m, quantity.field="Ms1.Area")





####################################
#now doing DEPC
####################################

df.spectr.H1 <- diann_load(here(input.dir,"Spectronaut_VISTAmodifications_DEPC(H)_Report.tsv"))



#reformat to get longer
df2<- df.spectr.H1 %>% select(!ends_with("EG.TotalQuantity (Settings)")) %>% 
	select(!contains("06302023")) %>%
	pivot_longer(cols=contains("DIA_iRT"), names_to="Run",values_to="MS2Quantity")

df2<-df2 %>% rename(`Protein.Group`="PG.ProteinGroups",
		#"Protein.Ids","Protein.Names",
		Genes="PG.Genes", 
		#"PG.MaxLFQ",
		#"Genes.MaxLFQ", "Modified.Sequence", 
		`Stripped.Sequence`="PEP.StrippedSequence", 
		`Precursor.Id`="EG.PrecursorId"
		#,"Precursor.Normalised", 
		#"Ms1.Area" 
		) %>% mutate(`Modified.Sequence`=`Precursor.Id`)

#parse out the condition/grouping info from the field of Run
df2$Run2<-sub(df2$Run, pattern=".raw.PEP.MS2Quantity",fixed=T,replacement="")
df2$Run2<-sub(df2$Run2, pattern="_[0-9]+$",replacement="") #<-get rid of date, trailing 
df2$Repeat<-substr(df2$Run2, nchar(df2$Run2), nchar(df2$Run2)) #<- get repeat number
df2$Run2<-sub(df2$Run2,pattern="\\-[1-3]{1}$", replacement="") #<- get rid of repeat number

df2$Run2<-sub(df2$Run2,pattern="^\\[[0-9]+\\] DIA_iRT_", replacement="") #<- get rid of leading "DIA_iRT"



#new remove other modifications (not DEPC)

#get rid of Unimod:4 which we don't want to take care of in this round
df2$Modified.Sequence <- sub(pattern="^_",
			x=df2$Modified.Sequence,replacement="") %>% 
			sub(pattern="_\\.[0-9]*$",
				replacement="") 
df2$Modified.Sequence <- 
		gsub(pattern="\\[(Carbamidomethyl|Oxidation) \\([CM]+\\)\\]",
			x=df2$Modified.Sequence,"")

##
# first we need to go through the modified sequences to get
# all possible modifications
mseq<-str_extract_all(pattern="\\[[a-zA-Z0-9]+\\s*\\([A-Z]+\\)\\]",
			string=df2$Modified.Sequence )

mods<-unique(unlist(mseq))

#with only DEPC(H) we can run

m<-mods[1]
out.file<-paste0(out.file0, "_",m,"_spectronaut_H1.csv")
df2$Modified.Sequence2<-df2$Modified.Sequence
		#get rid of mod[1] "(DEPC)", since we want to study "Fumi"
		# don't need to do this now.
		#  Don't RUN
			df2$Modified.Sequence2<-gsub(pattern=m_remove,
						x=df2$Modified.Sequence,"", fixed=T)
		##################end of don't run section###########
#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]
runModStats(dt=df2, output.dir, out.file, m, quantity.field="MS2Quantity")
