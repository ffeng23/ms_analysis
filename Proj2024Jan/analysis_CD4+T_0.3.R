# R code to do analysis of MS data with CD8+ T cell in vitro
# do statistics
#  1/9/2024
#  		copied this file from ../Proj3/analysis_CD4+T_0.3.R
# in this module, we do analysis
# 1) statistical analysis at each time point to compare between
# 		groups
# 2) do heatmeap and line plot to see patterns
# 		added line plot to the pattern
#
#-------------old comments from previous file------ 
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
library(openxlsx)
library(rstatix)


data.dir<-"Proj2024Jan"
output.dir<-"Proj2024Jan"

source(here("functions.R"))


#read the data
df <- diann_load(here(data.dir,
	"CD8+T_VISTA-SNS101_12292023_Report.tsv"))

#####---can not do it this time, since it missing peptid id
#disable this for now
#df<-removeReplicates(df) 

#pivot the table back to make the flow as before.
df.long<- df %>% #dplyr::select(!starts_with("HeLa")) %>%
	pivot_longer(cols=contains("CD8+T"),
	names_to="Run", values_to = "PG.MaxLFQ"
	)
#cleanup the record and keep for now only PG.MaxLFQ
# Run, 
df.clean <- df.long %>% dplyr::select(c("Run","PG.ProteinGroups",
		#"Protein.Ids", #"Protein.Names",
		"PG.Genes",	
		#"Precursor.Id",
	 "PG.MaxLFQ")) %>%
	mutate(Run2=sub(x=Run, pattern="_[0-9]+\\.raw\\.PG\\.Quantity$","")) %>%
	mutate(Run2=sub(x=Run2, pattern="^\\[[0-9]+\\] CD8\\+T_","")) %>%
	#mutate(Run2=sub(x=Run2, pattern="CD4+T_","", fixed=T)) %>%
	#mutate(Group="Reactome") %>%
	#mutate(Group=replace(Group, grepl(Run2, pattern="^Secret"),"Secrete")) %>%
	#mutate(Run2=sub(x=Run2, pattern="^Secret_","")) %>%
	mutate(Treatment=sub(x=Run2, pattern="^[0-9A-Za-z]+_","")) %>%
	mutate(Repeat=sub(x=Treatment, pattern="^(Contr|VISTA|VISTA-SNS-101)\\-","")) %>%
	
	mutate(Treatment=sub(x=Treatment, pattern="\\-[1-3]{1}$","")) %>%
	
	mutate(Run2=sub(x=Run2, pattern="\\-[1-3]{1}$","")) %>% 
	
	mutate(Time=sub(x=Run2, pattern="_(Contr|VISTA|VISTA-SNS-101)$","")) %>%
	
	mutate(Time=sub(x=Time,pattern="(hr|min)","")) %>%
	mutate(Time=replace(Time, grepl(Time,pattern="30"),"0.5"))%>%
	mutate(Time=as.numeric(Time))

df.clean <- df.clean %>% dplyr::rename(Protein.Group=PG.ProteinGroups,
		Genes=PG.Genes)
#now we need to split to have twoo groups for doing stats
#let's do average the repeat first (collapse the 3 repeats)
##NO!!! we can not do this, since there is only technical repeats
# split groups
df.clean<-df.clean %>% replace(is.na(.),0)

# do analysis for df.proteome
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
	mutate(y=log2(y+1))

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
	t_test(y ~ Treatment,detailed=T) %>% #View()
	adjust_pvalue(method="BH") %>%
  	add_significance("p")
)

saveRDS(df.stat, file=here(output.dir,"summary_t_test_CD8+T_fill_pivotData.Rds"))
#df.stat<-readRDS( file=here(output.dir,"summary_t_test_CD8+T_fill_pivotData.Rds"))

df.sumstat<- items.stat2 %>%
	group_by(Protein.Group, Genes,  Time, Treatment ) %>%
	get_summary_stats(y, type="common") 


saveRDS(df.sumstat, file=here(output.dir,"summary_stat_CD8+T_fill_pivotData.Rds"))
#df.sumstat<-readRDS(file=here(output.dir,"summary_stat_CD8+T_fill_pivotData.Rds"))

df.sumstat <- df.sumstat %>% dplyr::select(
	Protein.Group,Genes,Time,Treatment,mean
	)

df.stat.m1<-df.stat %>% left_join(df.sumstat, 
	by=c("Protein.Group","Genes", "Time", "group1" = "Treatment")) %>% 
	rename(mean1=mean)
df.stat.m2<-df.stat.m1 %>% left_join(df.sumstat, 
	by=c("Protein.Group","Genes", "Time", 
		"group2" = "Treatment")) %>% 
	rename(mean2=mean)

saveRDS(df.stat.m2, file=here(output.dir,"stat_CD8+T_fill_pivotData.Rds"))
#df.stat.m2<- readRDS(file=here(output.dir,"stat_CD8+T_fill_pivotData.Rds"))

write_csv(df.stat.m2, 
	file=here(output.dir,"stat_CD8+T_fill_pivotData.csv"))

#rewrite it into wider 
df.stat.m2.wider <- df.stat.m2 %>% 
	mutate(Time=paste0(Time,"hrs")) %>%
	rename(foldchange=estimate) %>%
	pivot_wider(id_cols=c("Protein.Group","Genes"),
			names_from=c(group1, group2,Time),
			values_from=c( #n1, n2, 
					foldchange,
					#statistic,df,p, 
					p.adj#, #p.signif,
					#p.adj.signif, mean1, mean2, 
				),
			names_vary="slowest",
			names_glue="{Time}_{group1}/{group2}_{.value}"
		) 
df.stat.m2.wider <- df.stat.m2.wider [,c(1,2,
		3,4,9,10, 15,16,21,22,27:28,
		5:6,11:12,17:18,23:24,29:30,
		7:8,13:14,19:20,25:26,31:32)]
write_csv(df.stat.m2.wider, 
	file=here(output.dir,"stat_CD8+T_fill_pivotData_wider_reformat.csv"))



#create summary statistics like means 
items.means<-	df.clean2 %>%
	#arrange(Stripped.Sequence, Time, Treatment)
	group_by(Protein.Group, Genes,  Time, Treatment ) %>%
	summarize(mean=mean(y), samplesize=n()) #%>% View("samplesize")

items.means.wide <- pivot_wider(items.means,
	id_cols=c("Protein.Group","Genes","Time"),
	names_from=c(Treatment), values_from=c(mean, samplesize))

items.all<- df.stat.m2 #df.stat %>% left_join(items.means.wide)

items.all.wide <- items.all %>% 
	mutate_at(c("group1", "group2"), as.factor) %>%
	pivot_wider(
		id_cols=c("Protein.Group","Genes"),
		names_from=c(Time,group1, group2), values_from=c(n1,n2, 
			statistic,p, p.adj,p.signif, p.adj.signif,
			mean1, mean2)

	)

#now let's do it by time for showing stats as output
#
times<-c(0.5, 2, 6, 12, 24)
wb <- createWorkbook("Feng")

for(ts in times)
{
	tmp<-items.all %>% 
		dplyr::filter(Time==ts) #%>% 
	addWorksheet(wb, sheetName=paste0(ts,"hrs"))
	writeData(wb, x=tmp, 
		sheet=paste0(ts,"hrs"))

}

saveWorkbook(wb,file=here(output.dir, "stats.xlsx"),
	 overwrite = TRUE) 

#now let's work on generating data to do trend,
# we basically split the data into vista-control
#  and vista-sns101 - control
trend.vista<-items.means.wide %>% 
	dplyr::select(!starts_with("samplesize")) %>%
	dplyr::select(!contains("VISTA-SNS101")) %>% 
	mutate(fc= mean_VISTA - mean_Contr)

trend.vista.wide<-trend.vista %>% pivot_wider(
		id_cols=c(Protein.Group,Genes), 
		values_from=fc, names_from=Time
	) %>% mutate("0"=0) %>% relocate("0",.before="0.5")
trend.vista.wide$"0"<-#trend.vista.wide$"0"+
	rnorm(length(trend.vista.wide$"0"),0,0.000001)

#save 
saveRDS(trend.vista.wide, 
	file=here(output.dir,"Trend.vista_vs_ctrl.RDS"))
	
pdf(file=here(output.dir,"VISTAeffect_raw_hier.pdf"),
	width=7, height=7)
v_h<-pheatmap(trend.vista.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE,# cluster_row=F
	labels_column=paste0(names(trend.vista.wide[,-c(1,2)]),"hrs"),
	show_rownames=F
	)
dev.off()
pdf(file=here(output.dir,"VISTAeffect_kmean.pdf"),
	width=7, height=4)
v_k<-pheatmap(trend.vista.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE, kmeans_k=18#, 
	#filename=here(output.dir,"VISTAeffect_kmean.pdf")
	)
dev.off()
trend.vista.wide$clusters=v_k$kmeans$cluster
clusters.vista<-trend.vista.wide[
		order(trend.vista.wide$clusters),
		c("Protein.Group","Genes","clusters")]
write_csv(clusters.vista,
	file=here(output.dir,"kmean_cluster_vista_k18.csv"))

#also draw line figures
trend.vista.long<- trend.vista.wide %>%
	pivot_longer(cols=c("0","0.5","2","6","12","24"),
		names_to="time_point",
		values_to="fc_over_zero"
		)

trend.vista.long$clusters<-factor(trend.vista.long$clusters,
		levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
trend.vista.long$time_point<-factor(trend.vista.long$time_point,
		levels=c("0","0.5","2","6","12","24"))

trend.vista.long.summary<-trend.vista.long %>%
	group_by(clusters, time_point) %>%
	summarize(fc_mean=mean(fc_over_zero),
		fc_median=median(fc_over_zero))

pdf(file=here(output.dir,"VISTAeffect_kmean_lines_k18.pdf"),
	width=10, height=6)
ggplot(		
	)+
	geom_line(data=trend.vista.long,aes(x=time_point, y=fc_over_zero, group=Protein.Group))+ xlab("Time Point (hr)")+ylab("FC over t=0hr")+
	geom_hline(yintercept=0,colour="grey",linetype=2)+
	geom_line(data=trend.vista.long.summary, 
		aes(x=time_point, y=fc_mean, group=clusters),
			color="red",linewidth=1.1, linetype=2)+
	facet_wrap(.~clusters)
dev.off()	


trend.vista_sns101<-items.means.wide %>% 
	dplyr::select(!starts_with("samplesize")) %>%
	dplyr::select(-`mean_VISTA`) %>% 
	mutate(fc= `mean_VISTA-SNS-101` - mean_Contr)

trend.vista_sns101.wide<-trend.vista_sns101 %>% pivot_wider(
		id_cols=c(Protein.Group,Genes), 
		values_from=fc, names_from=Time
	) %>% mutate("0"=0) %>% relocate("0",.before="0.5")
trend.vista_sns101.wide$"0"<-#trend.vista.wide$"0"+
	rnorm(length(trend.vista.wide$"0"),0,0.000001)

saveRDS(trend.vista_sns101.wide, 
	file=here(output.dir,"Trend.vistaSns101_vs_ctrl.RDS"))
pdf(file=here(output.dir,"VISTA_SNS101_effect_hier.pdf"),
	width=7, height=7)
pheatmap(trend.vista_sns101.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE, #cluster_row=F
	labels_column=paste0(names(trend.vista_sns101.wide[,-c(1,2)]),"hrs"),
	show_rownames=F
	)
dev.off()
pdf(file=here(output.dir,"VISTA_SNS101_effect_kmean.pdf"),
	width=7, height=4)
vs_k<-pheatmap(trend.vista_sns101.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE, kmeans_k=18
	)
dev.off()

trend.vista_sns101.wide$clusters=vs_k$kmeans$cluster
clusters.vista_sns101<-trend.vista_sns101.wide[
		order(trend.vista_sns101.wide$clusters),
		c("Protein.Group","Genes","clusters")]
write_csv(clusters.vista_sns101,
	file=here(output.dir,"kmean_cluster_vistasns101_k18.csv"))

#also draw line figures
trend.vista_sns101.long<- trend.vista_sns101.wide %>%
	pivot_longer(cols=c("0","0.5","2","6","12","24"),
		names_to="time_point",
		values_to="fc_over_zero"
		)

trend.vista_sns101.long$clusters<-factor(trend.vista_sns101.long$clusters,
		levels=c(1,2,3,4,5,6,7,8,9,10,11:18))
trend.vista_sns101.long$time_point<-factor(trend.vista_sns101.long$time_point,
		levels=c("0","0.5","2","6","12","24"))

trend.vista_sns101.long.summary<-trend.vista_sns101.long %>%
	group_by(clusters, time_point) %>%
	summarize(fc_mean=mean(fc_over_zero),
		fc_median=median(fc_over_zero))

pdf(file=here(output.dir,"VISTA_sns101effect_kmean_lines_k18.pdf"),
	width=10, height=6)
ggplot(		
	)+
	geom_line(data=trend.vista_sns101.long,aes(x=time_point, y=fc_over_zero, group=Protein.Group))+ xlab("Time Point (hr)")+ylab("FC over t=0hr")+
	geom_hline(yintercept=0,colour="grey",linetype=2)+
	geom_line(data=trend.vista_sns101.long.summary, 
		aes(x=time_point, y=fc_mean, group=clusters),
			color="red",linewidth=1.1, linetype=2)+
	facet_wrap(.~clusters)
dev.off()	



		#### sns101 vs vista
trend.vista_vs_sns101<-items.means.wide %>% 
	dplyr::select(!starts_with("samplesize")) %>%
	dplyr::select(-`mean_Contr`) %>% 
	mutate(fc= `mean_VISTA-SNS-101` - mean_VISTA)

trend.vista_vs_sns101.wide<-trend.vista_vs_sns101 %>% 
	pivot_wider(
		id_cols=c(Protein.Group,Genes), 
		values_from=fc, names_from=Time
	) %>% mutate("0"=0) %>% relocate("0",.before="0.5")

trend.vista_vs_sns101.wide$"0"<-#trend.vista.wide$"0"+
	rnorm(length(trend.vista_vs_sns101.wide$"0"),0,0.000001)

saveRDS(trend.vista_vs_sns101.wide, 
	file=here(output.dir,"Trend.vista_vs_vistaSNS101.RDS"))

pdf(file=here(output.dir,"VISTA_vs_SNS101_effect_hier.pdf"),
	width=7, height=7)
pheatmap(trend.vista_vs_sns101.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE, #cluster_row=F
	labels_column=paste0(names(trend.vista_vs_sns101.wide[,-c(1,2)]),"hrs"),
	show_rownames=F
	)
dev.off()
pdf(file=here(output.dir,"VISTA_vs_SNS101_effect_kmean.pdf"),
	width=7, height=4)
vs_k<-pheatmap(trend.vista_vs_sns101.wide[,-c(1,2)], scale="none",
	cluster_cols=FALSE, kmeans_k=15
	)
dev.off()

trend.vista_vs_sns101.wide$clusters=vs_k$kmeans$cluster
clusters.vista_vs_sns101<-trend.vista_vs_sns101.wide[
		order(trend.vista_vs_sns101.wide$clusters),
		c("Protein.Group","Genes","clusters")]
write_csv(clusters.vista_vs_sns101,
	file=here(output.dir,"kmean_cluster_vistasns101_vs_vista_k15.csv"))

#also draw line figures
trend.vista_vs_sns101.long<- trend.vista_vs_sns101.wide %>%
	pivot_longer(cols=c("0","0.5","2","6","12","24"),
		names_to="time_point",
		values_to="fc_over_zero"
		)

trend.vista_vs_sns101.long$clusters<-factor(trend.vista_vs_sns101.long$clusters,
		levels=c(1,2,3,4,5,6,7,8,9,10,11:15))
trend.vista_vs_sns101.long$time_point<-factor(trend.vista_vs_sns101.long$time_point,
		levels=c("0","0.5","2","6","12","24"))

trend.vista_vs_sns101.long.summary<-trend.vista_vs_sns101.long %>%
	group_by(clusters, time_point) %>%
	summarize(fc_mean=mean(fc_over_zero),
		fc_median=median(fc_over_zero))

pdf(file=here(output.dir,"VISTA_vs_sns101effect_kmean_lines_k15.pdf"),
	width=10, height=6)
ggplot(		
	)+
	geom_line(data=trend.vista_vs_sns101.long,aes(x=time_point, y=fc_over_zero, group=Protein.Group))+ xlab("Time Point (hr)")+ylab("FC over t=0hr")+
	geom_hline(yintercept=0,colour="grey",linetype=2)+
	geom_line(data=trend.vista_vs_sns101.long.summary, 
		aes(x=time_point, y=fc_mean, group=clusters),
			color="red",linewidth=1.1, linetype=2)+
	facet_wrap(.~clusters)
dev.off()	





	