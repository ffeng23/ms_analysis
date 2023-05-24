# R code to do analysis of MS data with T cell in vitro
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
#read the data
df <- diann_load(here("Data","DIA-NN.csv"))
df<-removeReplicates(df)


# replace na with 0
df[is.na(df)]<-0
df<-df %>% mutate_if( is.numeric, ~log2(.+1) )



#clean up the column names
group_info<- df %>%
	select(Protein.Group, Protein.Ids, Genes)
group_info[duplicated(group_info$Protein.Ids),"Protein.Ids"]<-
	paste0(group_info[duplicated(group_info$Protein.Ids),"Protein.Ids"],"_2")
group_info<-as.data.frame(group_info)
rownames(group_info)<-group_info$Protein.Ids
mat<-df %>% 
	select(-Protein.Group, -Protein.Ids, -Genes)

#clean up 
nms<-names(mat)
nms<-sub(x=nms, pattern="^CD4\\+T__iRT1000ng\\-DIA\\-",
		replace="Ctrl",fixed=F)

nms<-sub(x=nms, pattern="CD4+T_iRT_1000ng-DIA-",
		replace="",fixed=T)
nms<-sub(x=nms, pattern="_[0-9]{8}$",
		replacement="")
nms<-gsub(x=nms, pattern="_",replacement="-")
nms<-sub(x=nms, pattern="SNS-101",
	replacement="S101")
names(mat)<-nms

#now summarize the technical repeat into one
i<-1
mat2<-data.frame()
for(i in 1:9)
{
	temp<-apply(mat[,(i-1)*3+c(1:3)],MARGIN=1,mean)
	if(i==1){
		mat2=as.data.frame(temp)
	} else {
		temp=as.data.frame(temp)
		mat2<-cbind(mat2,temp)
	}
}
nms<-sub(x=nms,pattern="\\-[0-3]{1}$",replacement="")
names(mat2)<-unique(nms)
groups<-data.frame(id=c(1:27),samples=nms)#,group)
groups$group<-sub(groups$samples,pattern="[123]{1}$",replacement="")
groups$group<-sub(groups$group, pattern="\\-$",replacement="", fixed=F)

#now start doing the sample annotation
fullname<-data.frame(FullRunName=names(df)[-c(1:3)])
groups<-cbind(fullname, groups)
groups$MS_batch<-rep(rep(c(1:3),rep(3,3)),3)
groups$sample_num<-sub(x=groups$samples, pattern="[0-9A-Za-z]+\\-", 
	replacement="")

groups$sample_num<-sub(x=groups$sample_num, pattern="[a-zA-Z]+", 
	replacement="")

groups$sample_num<-sub(x=groups$sample_num, pattern="\\-", 
	replacement="")
groups$sample_num<-as.numeric(groups$sample_num)
groups$group_num<-rep(c(1:3),rep(9,3))
groups$repeat_num<-rep(c(1:3),9)
groups$order<-(groups$MS_batch-1)*9+((groups$group_num-1)*3+groups$repeat_num)

groups$MS_batch<-as.factor(groups$MS_batch)


#now let's draw pca
#doiing PCA without supplementary

res<-prcomp(t(mat),scale.=T, center=T)
fig_eig.res<-fviz_eig(res,addlabels = TRUE, ylim = c(0, 30))

pind_g<-fviz_pca_ind(res, axes=c(1,2),
             habillage = as.factor(groups$group), # color by groups 
			 geom=c("point","text"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = F, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
			 mean.point=FALSE,
             repel = TRUE # Avoid text overlapping
             )+ guides(shape=guide_legend(title="Group"),
             	color=guide_legend(title="Group")
             )
pind<-fviz_pca_ind(res, axes=c(1,2),
             habillage = as.factor(groups$MS_batch), # color by groups 
			 geom=c("point","text"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = F, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
			 mean.point=FALSE,
             repel = TRUE # Avoid text overlapping
             )+ guides(shape=guide_legend(title="Batch"),
             	color=guide_legend(title="Batch")
             )

#clearly, we need to do normalization. center 
# we need to do box plot for each sample
pdf(file=here("Proj2","normalization.pdf"))      
boxplot(mat)                       
dev.off()

pdf(file=here("Proj2","pca_raw.pdf"),width=15, height=6)      
ggarrange(pind,pind_g, ncol=2)                       
dev.off()

	###################################
	#now need to start using probatch for correction
	####################################

library(proBatch)

#define the data input, we will use df to start with and turn it into long form
df[2396,"Protein.Ids"]<-paste0(df[2396,"Protein.Ids"],"-b")

names(df)[2]<-'peptide_group_label'
df_long<- df %>% pivot_longer(cols=ends_with("12212022"),
		names_to="FullRunName", values_to="Intensity")

feature_id_col = 'peptide_group_label'
measure_col = 'Intensity'
sample_id_col = 'FullRunName'
essential_columns = c(feature_id_col, measure_col, sample_id_col)

technical_factors = c('MS_batch')
biological_factors = c('group')
biospecimen_id_col = 'samples'

batch_col='MS_batch'

## create feature annotation
#feature info is group_info 


# follow the code in tutorial/vignette
example_proteome = df_long %>% select(one_of(essential_columns))
gc()

example_matrix <- 
 long_to_matrix(example_proteome,
	feature_id_col = 'peptide_group_label',
	measure_col = 'Intensity',
	sample_id_col = 'FullRunName')

plot_sample_mean(example_matrix, groups, order_col = 'order',
batch_col = batch_col, color_by_batch = TRUE, ylimits = c(15, 18.5),
#color_scheme = color_list[[batch_col]]
)

batch_col = 'MS_batch'
bp<-plot_boxplot(example_proteome, groups,
batch_col = batch_col#, color_scheme = color_list[[batch_col]]
)
pdf(file=here("Proj2","boxplot_raw.pdf"),
		width=7, height=3)
bp 
dev.off()
quantile_normalized_matrix = normalize_data_dm(example_matrix,
normalize_func = 'quantile')

plot_sample_mean(quantile_normalized_matrix, groups,
color_by_batch = TRUE, ylimits = c(15, 18),
#color_scheme = color_list[[batch_col]]
)

selected_annotations <- c('MS_batch',
'group')

#Plot clustering between samples
plot_hierarchical_clustering(quantile_normalized_matrix,
	sample_annotation = groups,
	#color_list = color_list,
	factors_to_plot = selected_annotations,
	distance = 'euclidean', agglomeration = 'complete',
	label_samples = FALSE)

pca1 = plot_PCA(quantile_normalized_matrix, groups, color_by = 'MS_batch',
plot_title = 'MS batch'#, color_scheme = color_list[['MS_batch']]
)+scale_size_manual(values = 15) 

pca2 = plot_PCA(quantile_normalized_matrix, groups, color_by = 'group',
plot_title = 'group'#, color_scheme = color_list[['MS_batch']]
)

ggarrange(pca1, pca2, ncol = 2, nrow = 1)

pvca1<-plot_PVCA(quantile_normalized_matrix, groups,
technical_factors = technical_factors,
biological_factors = biological_factors)


quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix)
loess_fit_70 <- adjust_batch_trend_df(quantile_normalized_long, 
	groups,
span = 1)

plot_with_fitting_curve(feature_name = 'A0A024RBG1',
fit_df = loess_fit_70, fit_value_col = 'fit',
df_long = quantile_normalized_long,
sample_annotation = groups, color_by_batch = TRUE, 
#color_scheme = color_list[[batch_col]],
plot_title = 'Span = 60%')


###correct with feature level 
#peptide_median_df <- center_feature_batch_medians_df(loess_fit_70, groups)

#peptide_median_df <- center_feature_batch_medians_df(quantile_normalized_long, groups)

peptide_median_df<-correct_with_ComBat_df(quantile_normalized_long, groups)
#peptide_median_df<-correct_with_ComBat_df(loess_fit_70, groups)

batch_corrected_matrix = long_to_matrix(peptide_median_df)

pca.corrected = plot_PCA(batch_corrected_matrix, groups, color_by = 'MS_batch',
plot_title = 'MS batch'#, color_scheme = color_list[['MS_batch']]
)

x<-apply(batch_corrected_matrix,1, FUN=function(x){length(unique(x))==1})
y<-which(x)
d<-dim(batch_corrected_matrix)[2]
batch_corrected_matrix[y,]<-batch_corrected_matrix[y,]+
	matrix(rnorm(d*length(y),0,0.01),nrow=length(y),
		ncol=d)
batch_corrected_matrix2<-batch_corrected_matrix
colnames(batch_corrected_matrix2)<-
	paste0(groups$group,'-',groups$MS_batch,'-',groups$repeat_num)
res2<-prcomp(t(batch_corrected_matrix2),scale.=T, center=T)

pind2<-fviz_pca_ind(res2, axes=c(1,2),
             habillage = as.factor(groups$MS_batch), # color by groups 
			 geom=c("point","text"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = F, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
			 mean.point=FALSE,
             repel = TRUE # Avoid text overlapping
             )+ guides(shape=guide_legend(title="Batch"),
             	color=guide_legend(title="Batch")
             )
	
pind2_g<-fviz_pca_ind(res2, axes=c(1,2),
             habillage = as.factor(groups$group), # color by groups 
			 geom=c("point","text"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = F, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
			 mean.point=FALSE,
             repel = TRUE # Avoid text overlapping
             )
			

pdf(file=here("Proj2","pca_corrected.pdf"),width=15, height=6)
ggarrange(pind2, pind2_g, ncol=2)
dev.off()

pca.corrected2 = plot_PCA(batch_corrected_matrix, groups, color_by = 'group',
plot_title = 'group'#, color_scheme = color_list[['MS_batch']]
)

ggarrange(pca.corrected, pca.corrected2, ncol = 2, nrow = 1)

pvca2<-plot_PVCA(batch_corrected_matrix, groups,
technical_factors = technical_factors,
biological_factors = biological_factors)
pdf(file=here("Proj2","PVCA_both.pdf"),width=8,height=7)
ggarrange(pvca1, pvca2, nrow=2)
dev.off()

#now let's do heatmap

plot_heatmap_diagnostic(batch_corrected_matrix2, groups,
factors_to_plot = selected_annotations,
cluster_cols = T,cluster_rows=T,
#color_list = color_list,
show_rownames = F, show_colnames = T)

library(pheatmap)

heat.all<-pheatmap(batch_corrected_matrix2,
		scale="row", cluster_rows=T,
		cluster_cols=F, 
		show_rownames=F)

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
library(rstatix)
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



#now start doing the stats. on all samples including technical repeats.
#
peptide_median_df.group.vi2 <- peptide_median_df.group
peptide_median_df.group.vi2$Intensity<-peptide_median_df.group.vi2$Intensity +
	rnorm(length(peptide_median_df.group.vi2$Intensity),0,0.01)
peptide_median_df.group.vi2 <- peptide_median_df.group.vi2 %>%	
	dplyr::filter(group!="S101-VISTA") %>% 
	group_by(peptide_group_label) %>% 
	t_test(Intensity ~ group) %>%
	adjust_pvalue() %>%
  	add_significance("p.adj")
