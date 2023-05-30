# R code to draw heatmaps
# using the previously saved data
#========================================
# updated 5/25/2023. rerun to use the removed duplicated data
# and redraw the heatmaps
#

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggpubr)
library(pheatmap)
library(readxl)

#define the data source
#   updated source to use the remove duplicated data
output.dir<-"Proj2/removeReps"
# to go back to original data please use below
#output.dir <-"Proj2"

#read data
df<-read_tsv(file=here(output.dir,"VistaCtrl_expression_BatchCorrected.tsv"))
df<- df %>% select(peptide_group_label,Ctrl1, Ctrl2, Ctrl3, "VISTA-1",
		"VISTA-2","VISTA-3", Genes)

gpath <- read_excel(here("Proj2","Gene symbols for heatmaps.xlsx"))
gpath<-gpath[,c(2,4,5)]
gpath<-as.data.frame(gpath)
#get the genes out
pathNames<-names(gpath)
names(gpath)<-c("tolls","TCRs","APCs")

#do first toll like receptors
tolls<-gpath$tolls
tolls<-tolls[1:23]

df.heat<- df %>% filter(is.element(Genes,tolls))
df.heat<-as.data.frame(df.heat)

df.heat<-df.heat[!duplicated(df.heat$Genes),]#<-
	#paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
#df.heat[duplicated(df.heat$Genes),"Genes"]<-
#	paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
rownames(df.heat)<-df.heat$Genes

df.heat.matrix<-df.heat %>% select(-Genes,-peptide_group_label)
#ann<-df.heat %>%

heat.toll<-pheatmap(df.heat.matrix,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=T)

##now let's do the second one.

tolls<-gpath$TCRs
tolls<-tolls[1:41]

df.heat<- df %>% filter(is.element(Genes,tolls))
df.heat<-as.data.frame(df.heat)

df.heat<-df.heat[!duplicated(df.heat$Genes),]#<-
	#paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
#df.heat[duplicated(df.heat$Genes),"Genes"]<-
#	paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
rownames(df.heat)<-df.heat$Genes

df.heat.matrix<-df.heat %>% select(-Genes,-peptide_group_label)
#ann<-df.heat %>%

heat.TCRs<-pheatmap(df.heat.matrix,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=T)

## do the third one
tolls<-gpath$APCs
tolls<-tolls[1:29]

df.heat<- df %>% filter(is.element(Genes,tolls))
df.heat<-as.data.frame(df.heat)

df.heat<-df.heat[!duplicated(df.heat$Genes),]#<-
	#paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
#df.heat[duplicated(df.heat$Genes),"Genes"]<-
#	paste0(df.heat[duplicated(df.heat$Genes),"Genes"],".")
rownames(df.heat)<-df.heat$Genes

df.heat.matrix<-df.heat %>% select(-Genes,-peptide_group_label)
#ann<-df.heat %>%

heat.APCs<-pheatmap(df.heat.matrix,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=T)

pdf(file=here(output.dir,"pathAPCs_heat.pdf"),
		width=4,height=4)
heat.APCs
dev.off()
pdf(file=here(output.dir,"pathTCRs_heat.pdf"),
		width=4,height=4.5)
heat.TCRs
dev.off()
pdf(file=here(output.dir,"pathTolls_heat.pdf"),
		width=4,height=3.0)
heat.toll
dev.off()








#########left over from previous one.

#summary to get mean of 3 reps of each sample using long format
peptide_median_df.group<-cbind(peptide_median_df,
		groups[peptide_median_df$FullRunName,-1])

peptide_df.sample<- peptide_median_df.group %>%
	select(samples, group, Intensity, peptide_group_label) %>%
	group_by(samples,group,peptide_group_label) %>% 
		summarize(mIntensity=mean(Intensity))

#now turn it back to matrix to do plotting.
peptide_df.sample.matrix <- as.data.frame(peptide_df.sample) %>%
	select(-group) %>%
	pivot_wider(
	names_from=samples,
	values_from=mIntensity
	)

peptide_df.sample.vi2<-peptide_df.sample
peptide_df.sample.vi2$mIntensity<-peptide_df.sample.vi2$mIntensity +
	rnorm(length(peptide_df.sample.vi2$mIntensity),0,0.01)
peptide_df.sample.vi2<-peptide_df.sample.vi2 %>%
	dplyr::filter(group!="S101-VISTA") %>%
	group_by(peptide_group_label) %>% 
	t_test(mIntensity ~ group) %>%
	#adjust_pvalue(method="BH") %>%
  	add_significance("p")
peptide_df.sample.vi2$Genes<-group_info[peptide_df.sample.vi2$peptide_group_label,"Genes"]
peptide_df.sample$Genes<-group_info[peptide_df.sample$peptide_group_label, "Genes"]
out<-peptide_df.sample.matrix
out$Genes<-group_info[out$peptide_group_label,"Genes"]

write_tsv(file=here("Proj2","VistaCtrl_stats.tsv"),
	peptide_df.sample.vi2)

write_tsv(file=here("Proj2","VistaCtrl_expression_BatchCorrected.tsv"),
	out)

peptide_df.sample.matrix2<-peptide_df.sample.matrix[,-1]

x<-apply(peptide_df.sample.matrix2,1, FUN=function(x){length(unique(x))==1})
y<-which(x)
d<-dim(peptide_df.sample.matrix2)[2]
peptide_df.sample.matrix2[y,]<-peptide_df.sample.matrix2[y,]+
	matrix(rnorm(d*length(y),0,0.01),nrow=length(y),
		ncol=d)


heat.sample<-pheatmap(peptide_df.sample.matrix2,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=F)


#now let's do only ctrl vs vista
peptide_df.mat.vista2<-peptide_df.sample.matrix2[,
	c("Ctrl1","Ctrl2", "Ctrl3", "VISTA-1","VISTA-2","VISTA-3")]

x<-apply(peptide_df.mat.vista2,1, FUN=function(x){length(unique(x))==1})
y<-which(x)
d<-dim(peptide_df.mat.vista2)[2]
peptide_df.mat.vista2[y,]<-peptide_df.mat.vista2[y,]+
	matrix(rnorm(d*length(y),0,0.01),nrow=length(y),
		ncol=d)
heat.vista2<-pheatmap(peptide_df.mat.vista2,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=F)
pdf(file=here("Proj2","Vista_vs_Ctrl_heat.pdf"),width=4,height=6)
heat.vista2
dev.off()

pdf(file=here("Proj2","all3_heat.pdf"),width=4,height=6)
heat.sample
dev.off()



#now start doing the stats.
peptide_median_df.group.vi2 <- peptide_median_df.group
peptide_median_df.group.vi2$Intensity<-peptide_median_df.group.vi2$Intensity +
	rnorm(length(peptide_median_df.group.vi2$Intensity),0,0.01)
peptide_median_df.group.vi2 <- peptide_median_df.group.vi2 %>%	
	dplyr::filter(group!="S101-VISTA") %>% 
	group_by(peptide_group_label) %>% 
	t_test(Intensity ~ group) %>%
	adjust_pvalue() %>%
  	add_significance("p.adj")
