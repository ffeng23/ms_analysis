# R code to do analysis of MS data with CD4+ T cell in vitro
# do statistics
#=================================
#	8/14/2023
# this really should be Kasumi...v0.3.R
# since it copies the Kasumi...v0.2.R
#	in this module, we do analayis of MS data with CD4+ T in vitro
#	1) statistical analysis at each time point with control vs. vista, control vs. sns101
#	2) then do heatmap to see patterns with all data and 
#
###################################
#    old comments     below           #
#----------------------------
#			8/4/2023
# #same data, but was combined according to
#	 PG group and peptides and counts, etc (Using unique values)
#	need to do removeDuplicate??
#
#------------------------------
#    7/30/2023
# see the document for description
#file:///home/feng/Feng/hg/ms_analysis/Proj3/Data/VISTAmodi & Kasumi3-SNS-101_timecourse for statistical analysis.pptx
#
library(here)
library(dplyr)
library(tidyr)
library(diann)
library(readr)
library(pheatmap)


data.dir<-"Proj3/Data"
output.dir<-"Proj3"

source("functions.R")
#read the data
df <- diann_load(here(data.dir,"CD4+T_VISTA-SNS101_Timecouse_07162023_report_PG.MaxLFQ_pivoted.csv"))
df<-removeReplicates(df) 

#pivot the table back to make the flow as before.
df.long<- df %>% dplyr::select(!starts_with("HeLa")) %>%
	pivot_longer(cols=starts_with("CD4+T"),
	names_to="Run", values_to = "PG.MaxLFQ"
	)
#cleanup the record and keep for now only PG.MaxLFQ
# Run, 
df.clean <- df.long %>% dplyr::select(c("Run","Protein.Group",
		"Protein.Ids", #"Protein.Names",
		"Genes",	
		#"Precursor.Id",
	 "PG.MaxLFQ")) %>%
	mutate(Run2=sub(x=Run, pattern="_[0-9]+$","")) %>%
	mutate(Run2=sub(x=Run2, pattern="^[a-zA-Z]+3_","")) %>%
	mutate(Run2=sub(x=Run2, pattern="CD4+T_","", fixed=T)) %>%
	#mutate(Group="Reactome") %>%
	#mutate(Group=replace(Group, grepl(Run2, pattern="^Secret"),"Secrete")) %>%
	#mutate(Run2=sub(x=Run2, pattern="^Secret_","")) %>%
	mutate(Treatment=sub(x=Run2, pattern="_[0-9A-Za-z\\-]+$","")) %>%
	mutate(Run2=sub(x=Run2, pattern="^(Control|VISTA|VISTA-SNS101)_","")) %>% 
	mutate(Time=sub(x=Run2, pattern="\\-[1-3]{1}$","")) %>%
	mutate(Repeat=sub(x=Run2, pattern="^[0-9]{1,2}(hr|min)\\-","")) %>%
	mutate(Time=sub(x=Time,pattern="(hr|min)","")) %>%
	mutate(Time=replace(Time, grepl(Time,pattern="30"),"0.5"))%>%
	mutate(Time=as.numeric(Time))
#now we need to split to have twoo groups for doing stats
#let's do average the repeat first (collapse the 3 repeats)
##NO!!! we can not do this, since there is only technical repeats
# split groups
df.clean<-df.clean %>% replace(is.na(.),0)

# do analysis for df.proteome
library(rstatix)
#now do t-test for each protein group
#
#df.proteome<-df.proteome %>% 
#	dplyr::filter(Time=="6hr"|Time=="24hr")
df.clean <-df.clean %>%
	mutate_at(c("Time","Treatment","Repeat"), as.factor)
#remeber to do log transformation with pseudo-count
df.clean2<- df.clean %>% 
	group_by( Protein.Group, Genes, Time, Treatment, Repeat,.drop=T) %>%
	summarize(y=mean(PG.MaxLFQ)) %>% 
	replace(is.na(.),0) %>%
	mutate(y=log(y+1))

#show the sample size and remove the one without repeats
items.stat<-	df.clean2 %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Protein.Group, Genes,  Time, Treatment ) %>%
	summarize(counts=n()) %>% #View("samplesize")
	dplyr::filter(counts>1)
# join/merge to have repeats (raw data to do stats)
items.stat2<- items.stat %>% left_join(df.clean2)

system.time(

df.stat<- items.stat2 %>%
	group_by(Protein.Group, Genes,  Time ) %>%
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")
)
saveRDS(df.stat, file=here(output.dir,"stat_CD4+T_fill_pivotData.Rds"))

#df.proteome.stat<-readRDS(file=here(output.dir,"stat_proteome_fill.Rds"))
write_csv(df.stat, 
	file=here(output.dir,"stat_CD4+T_fill_pivotData.csv"))





df.proteome.stat.wider<-pivot_wider(df.proteome.stat,id_cols=c("Protein.Group","Genes"),
	names_from=Time, values_from=c(n1,n2, statistic,p, p.adj,p.signif))
write_csv(df.proteome.stat.wider, 
	file=here(output.dir,"stat_proteome_fill_reformatted_pivotData.csv"))
#put the data into wider format

# do analysis for df.secreteome

#now do t-test for each protein group
#
#df.proteome<-df.proteome %>% 
#	dplyr::filter(Time=="6hr"|Time=="24hr")
df.secreteome <-df.secreteome %>%
	mutate_at(c("Time","Treatment","Repeat"), as.factor)
df.secreteome.stat<- df.secreteome %>% 
	group_by(Protein.Group, Genes, Time, Treatment, Repeat,.drop=FALSE) %>%
	summarize(y=mean(PG.MaxLFQ)) %>% 
	replace(is.na(.),0) %>%
	mutate(y=log(y+1))%>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by( Protein.Group, Genes,  Time ) %>%
	#summarize(counts=n())
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")

saveRDS(df.secreteome.stat, file=here(output.dir,"stat_secreteome_fill_pivotData.Rds"))
#df.secreteome.stat<-readRDS(file=here(output.dir,"stat_secreteome_fill.Rds"))


write_csv(df.secreteome.stat, 
	file=here(output.dir,"stat_secrete_fill_pivotData.csv"))

#reformat to wider 
df.secreteome.stat.wider<-pivot_wider(df.secreteome.stat,id_cols=c("Protein.Group","Genes"),
	names_from=Time, values_from=c(n1,n2, statistic,p, p.adj,p.signif))
write_csv(df.secreteome.stat.wider, 
	file=here(output.dir,"stat_secreteome_fill_reformatted_pivotData.csv"))



#now here we remove the one with too few repeats due 
# to undetected.
df.proteome.list<- df.proteome %>% 
	#group_by(Stripped.Sequence, Protein.Group, Time, Treatment, Repeat,.drop=FALSE) %>%
	#summarize(y=mean(PG.MaxLFQ)) %>% 
	#replace(is.na(.),0) %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Protein.Group,Genes, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% #View("f") 
	dplyr::filter(counts<2)

df.proteome2<-df.proteome %>%
	dplyr::filter(!is.element(Protein.Group,
			unique(df.proteome.list$Protein.Group )))

df.show<-df.proteome2 %>% 
	group_by(Protein.Group,Genes,  Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% View("f") 
	#dplyr::filter(counts<3)

df.proteome.stat2<- df.proteome2 %>% 
	mutate(y=log(PG.MaxLFQ+1)) %>%
	group_by(Protein.Group, Genes,  Time)%>%
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")

# df.proteome2[df.proteome2$Stripped.Sequence=="AAAAAAGAGPEMVR",]
saveRDS(df.proteome.stat2, file=here(output.dir,"stat_proteome_conserved_pivotData.Rds"))
#df.proteome.stat2<-readRDS(file=here(output.dir,"stat_proteome_conserved.Rds"))
write_csv(df.proteome.stat2, 
	file=here(output.dir,"stat_proteome_conserved_pivotData.csv"))

#reformat to wider 
df.proteome.stat2.wider<-pivot_wider(df.proteome.stat2,id_cols=c("Protein.Group","Genes"),
	names_from=Time, values_from=c(n1,n2, statistic,p, p.adj,p.signif))
write_csv(df.proteome.stat2.wider, 
	file=here(output.dir,"stat_proteome_conserved_reformatted_pivotData.csv"))




#now here we remove the one with too few repeats due 
# to undetected.
df.secreteome.list<- df.secreteome %>% 
	#group_by(Stripped.Sequence, Protein.Group, Time, Treatment, Repeat,.drop=FALSE) %>%
	#summarize(y=mean(PG.MaxLFQ)) %>% 
	#replace(is.na(.),0) %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Protein.Group, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% #View("f") 
	dplyr::filter(counts<2)

df.secreteome2<-df.secreteome %>%
	dplyr::filter(!is.element(Protein.Group,
			unique(df.secreteome.list$Protein.Group )))

df.show<-df.secreteome2 %>% 
	group_by(Protein.Group,Genes,  Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% View("f") 
	#dplyr::filter(counts<3)

df.secreteome.stat2<- df.secreteome2 %>% 
	mutate(y=log(PG.MaxLFQ+1)) %>%
	group_by(Protein.Group, Genes,  Time)%>%
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")


saveRDS(df.secreteome.stat2, file=here(output.dir,"stat_secrete_conserved_pivotData.Rds"))
#df.secreteome.stat2<-readRDS(file=here(output.dir,"stat_secrete_conserved.Rds"))
write_csv(df.secreteome.stat2, 
	file=here(output.dir,"stat_secrete_conserved_pivotData.csv"))

#reformat to wider 
df.secreteome.stat2.wider<-pivot_wider(df.secreteome.stat2,id_cols=c("Protein.Group","Genes"),
	names_from=Time, values_from=c(n1,n2, statistic,p, p.adj,p.signif))
write_csv(df.secreteome.stat2.wider, 
	file=here(output.dir,"stat_secreteome_conserved_reformatted_pivotData.csv"))

