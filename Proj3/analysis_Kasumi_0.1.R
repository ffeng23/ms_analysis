# R code to do analysis of MS data with kasumi T cell in vitro
# do statistics
# see the document for description
#file:///home/feng/Feng/hg/ms_analysis/Proj3/Data/VISTAmodi & Kasumi3-SNS-101_timecourse for statistical analysis.pptx
#
library(here)
library(dplyr)
library(tidyr)
library(diann)
library(readr)


data.dir<-"Proj3/Data"
output.dir<-"Proj3"


#read the data
df <- diann_load(here(data.dir,"Kasumi3_SNS-101_Timecourse_05272023.tsv"))
#df<-removeReplicates(df) 


#cleanup the record and keep for now only PG.MaxLFQ
# Run, 
df.clean <- df %>% select(c("Run","Protein.Group",
		"Protein.Ids", "Protein.Names","Genes",	
		"Precursor.Id",
	 "PG.MaxLFQ")) %>%
	mutate(Run2=sub(x=Run, pattern="_[0-9]+$","")) %>%
	mutate(Run2=sub(x=Run2, pattern="^[a-zA-Z]+3_","")) %>%
	mutate(Group="Reactome") %>%
	mutate(Group=replace(Group, grepl(Run2, pattern="^Secret"),"Secrete")) %>%
	mutate(Run2=sub(x=Run2, pattern="^Secret_","")) %>%
	mutate(Treatment=sub(x=Run2, pattern="_[0-9A-Za-z\\-]+$","")) %>%
	mutate(Run2=sub(x=Run2, pattern="^[0-9a-zA-Z]+_","")) %>% 
	mutate(Time=sub(x=Run2, pattern="\\-[1-3]{1}$","")) %>%
	mutate(Repeat=sub(x=Run2, pattern="^[0-9]{1,2}hr\\-",""))

#now we need to split to have twoo groups for doing stats
#let's do average the repeat first (collapse the 3 repeats)
##NO!!! we can not do this, since there is only technical repeats
# split groups
df.proteome<-df.clean %>% 
	filter(Group=="Reactome")

df.secreteome<-df.clean %>% 
	filter(Group=="Secrete")

# do analysis for df.proteome
library(rstatix)
#now do t-test for each protein group
#
df.proteome<-df.proteome %>% 
	dplyr::filter(Time=="6hr"|Time=="24hr")
df.proteome <-df.proteome %>%
	mutate_at(c("Time","Treatment","Repeat"), as.factor)
#remeber to do log transformation with pseudo-count
system.time(
df.proteome.stat<- df.proteome %>% 
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names,Time, Treatment, Repeat,.drop=FALSE) %>%
	summarize(y=mean(PG.MaxLFQ)) %>% 
	replace(is.na(.),0) %>%
	mutate(y=log(y+1))%>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names, Time ) %>%
	#summarize(counts=n())

	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")
)
saveRDS(df.proteome.stat, file=here(output.dir,"stat_proteome_fill.Rds"))

write_csv(df.proteome.stat, 
	file=here(output.dir,"stat_proteome_fill.csv"))
# do analysis for df.secreteome

#now do t-test for each protein group
#
#df.proteome<-df.proteome %>% 
#	dplyr::filter(Time=="6hr"|Time=="24hr")
df.secreteome <-df.secreteome %>%
	mutate_at(c("Time","Treatment","Repeat"), as.factor)
df.secreteome.stat<- df.secreteome %>% 
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names,Time, Treatment, Repeat,.drop=FALSE) %>%
	summarize(y=mean(PG.MaxLFQ)) %>% 
	replace(is.na(.),0) %>%
	mutate(y=log(y+1))%>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names, Time ) %>%
	#summarize(counts=n())
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")

saveRDS(df.secreteome.stat, file=here(output.dir,"stat_secreteome_fill.Rds"))
write_csv(df.secreteome.stat, 
	file=here(output.dir,"stat_secrete_fill.csv"))

#now here we remove the one with too few repeats due 
# to undetected.
df.proteome.list<- df.proteome %>% 
	#group_by(Stripped.Sequence, Protein.Group, Time, Treatment, Repeat,.drop=FALSE) %>%
	#summarize(y=mean(PG.MaxLFQ)) %>% 
	#replace(is.na(.),0) %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Precursor.Id, Protein.Group, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% #View("f") 
	dplyr::filter(counts<2)

df.proteome2<-df.proteome %>%
	dplyr::filter(!is.element(Precursor.Id,
			unique(df.proteome.list$Precursor.Id )))

df.show<-df.proteome2 %>% 
	group_by(Precursor.Id, Protein.Group,Genes, Protein.Names, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% View("f") 
	#dplyr::filter(counts<3)

df.proteome.stat2<- df.proteome2 %>% 
	mutate(y=log(PG.MaxLFQ+1)) %>%
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names, Time)%>%
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")

# df.proteome2[df.proteome2$Stripped.Sequence=="AAAAAAGAGPEMVR",]
saveRDS(df.proteome.stat2, file=here(output.dir,"stat_proteome_conserved.Rds"))

write_csv(df.proteome.stat2, 
	file=here(output.dir,"stat_proteome_conserved.csv"))

#now here we remove the one with too few repeats due 
# to undetected.
df.secreteome.list<- df.secreteome %>% 
	#group_by(Stripped.Sequence, Protein.Group, Time, Treatment, Repeat,.drop=FALSE) %>%
	#summarize(y=mean(PG.MaxLFQ)) %>% 
	#replace(is.na(.),0) %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Precursor.Id, Protein.Group, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% #View("f") 
	dplyr::filter(counts<2)

df.secreteome2<-df.secreteome %>%
	dplyr::filter(!is.element(Precursor.Id,
			unique(df.secreteome.list$Precursor.Id )))

df.show<-df.secreteome2 %>% 
	group_by(Precursor.Id, Protein.Group,Genes, Protein.Names, Time,Treatment,.drop=FALSE ) %>%
	summarize(counts=n()) %>% View("f") 
	#dplyr::filter(counts<3)

df.secreteome.stat2<- df.secreteome2 %>% 
	mutate(y=log(PG.MaxLFQ+1)) %>%
	group_by(Precursor.Id, Protein.Group, Genes, Protein.Names, Time)%>%
	t_test(y ~ Treatment) %>%
	adjust_pvalue(method="BH") %>%
  	add_significance("p")


saveRDS(df.secreteome.stat2, file=here(output.dir,"stat_secrete_conserved.Rds"))

write_csv(df.secreteome.stat2, 
	file=here(output.dir,"stat_secrete_conserved.csv"))
