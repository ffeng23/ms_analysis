# R code to do analysis of MS data with T cell in vitro
#
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

# output folder for new results after remove duplicates
# 
output.dir<-"Proj3"
input.dir<-"Proj3/Data"
out.file<-"VISTA_modifications_stats"

#read the data
df <- diann_load(here(input.dir,"VISTA_modifications_report - Copy.csv"))

df2<- df %>% dplyr::select("Run", "Protein.Group","Protein.Ids","Protein.Names","Genes", "PG.MaxLFQ",
		"Genes.MaxLFQ", "Modified.Sequence", "Stripped.Sequence", "Precursor.Id","Precursor.Normalised", "Ms1.Area" )
df2$Run2<-sub(df2$Run, pattern="_[0-9]+$",replacement="")
df2$Repeat<-substr(df2$Run2, nchar(df2$Run2), nchar(df2$Run2))
df2$Run2<-sub(df2$Run2,pattern="\\-[1-3]{1}$", replacement="")

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
m_remove<-mods[2]
m<-mods[1]
out.file<-paste0(out.file, "_",m,".csv")

#get rid of mod[1] "(DEPC)", since we want to study "Fumi"
df2$Modified.Sequence2<-gsub(pattern=m_remove,
			x=df2$Modified.Sequence,"", fixed=T)

#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]

sink(file=here(output.dir,out.file),append=F)
cat("Modification stat for ",m ,"\n")
sink()
#ss<-unique(df2$Stripped.Sequence)[1]
#ss<-"FKVATPYSLYVCPEGQNVTLTCR" 
#ss<-"FKVATPYSLYVCPEGQNVTLTCRLLGPVDK"
for(ss in unique(df2$Stripped.Sequence))
{
	cat("doing peptide sequence:",ss,".....\n")
	sink(file=here(output.dir,out.file),append=T)
	cat("'=================\n")
	cat("Peptide sequence ",ss ,"\n")
	sink()
	temp<-which(df2$Stripped.Sequence==ss)
	temp_df<-df2[temp,]
	#we need to figure out the unique different modification on different
	# amino acids
	m_indices<-gregexpr(text=temp_df$Modified.Sequence2, pattern="(",fixed=T)
	
	#now we need to figure out the stripped sequence index of modified aa
	m_indices_raw<-lapply(m_indices, FUN=
		function(x){
			y=x
			if(x[1]!=-1){
				y=(c(1:length(x))-1)*nchar(m)  # count for each position what is the raw index by subtracting the modified string
				y=x-y-1 # minus 1 because we are count at this point the starting position of "("
			}
			return (y)
		}
	 )

	m_indices_raw_unique<-unique(unlist(m_indices_raw))
	m_indices_raw_unique<-m_indices_raw_unique[-which(m_indices_raw_unique==-1)]	
	m_indices_raw_unique<-sort(m_indices_raw_unique)

	#for each we generate a stats 
	#m_raw=21
	for(m_raw in m_indices_raw_unique)
	{
		
		#get the string with modification *only*
		m_string<-paste0(substr(ss,1,m_raw), m, substr(ss, m_raw+1, nchar(ss)))
		m_string_start<-paste0(substr(ss,1,m_raw), m)
		cat("\t**doing modified sequence:", m_string_start,"\n")
		#get stats
		temp_df2<-temp_df %>% 
			mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T))%>%
			mutate_at(c("Run2","Run","Repeat"), as.factor) 
			
		stats <- temp_df2 %>% 
			#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
			mutate_at(c("Run2","Run","Repeat", "Modified"), as.factor) %>% 
			group_by(Run2, Repeat, Modified, .drop=F) %>% 
			summarize(Mod_Level=sum(Ms1.Area)
				)
		stats2 <- stats %>% 
			#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
			#mutate_at(c("Run2","Run","Repeat"), as.factor) %>% 
			group_by(Run2, Repeat) %>% 
			mutate(pct=Mod_Level/sum(Mod_Level)*100)  %>% ungroup

		stats3 <- stats2 %>% dplyr::filter(Modified==TRUE)%>%
			group_by (Run2) %>% 
			summarize(mean=mean(pct), std=sd(pct)) 

		sink(file=here(output.dir,out.file),append=T)
			cat("'***************\n")
			cat("\tModification at position ",m_raw ,"\n")
			cat("\tModified sequence:", m_string,"\n")
			write.csv(stats2)
			cat("\tnote: Mod_Level is the sum of Ms1.Area; pct is % of the total.\n")
			write.csv(stats3)
			cat("\tnote: mean is the mean % of 3 replicates; and sd is the standard deviation.\n")

		sink()
		#this last step is very important: remove the modifications 
		# this will make sure the patter in the following step 
		# can be matched.
		temp_df$Modified.Sequence2 <-sub(x=temp_df$Modified.Sequence2, 
				pattern=m_string_start,replacement=substr(ss,1,m_raw), fixed=T)
	}
}
