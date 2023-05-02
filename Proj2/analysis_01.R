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

#read the data
df <- diann_load(here("Data","DIA-NN.csv"))

# replace na with 0
df[is.na(df)]<-0
df<-df %>% mutate_if( is.numeric, ~log2(.+1) )

#clean up the column names
group_info<- df %>%
	select(Protein.Group, Protein.Ids, Genes)

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
#now let's draw pca
#doiing PCA without supplementary

res<-prcomp(t(mat),scale.=T, center=T)
fig_eig.res<-fviz_eig(res,addlabels = TRUE, ylim = c(0, 30))

pind<-fviz_pca_ind(res, axes=c(1,2),
             habillage = as.factor(groups$group), # color by groups 
			 geom=c("point","text"), #choice=c("group"),
             #palette = c("#00AFBB", "#E7B800","red","blue","pink","),
             addEllipses = F, #ellipse.type = "confidence", 
			 ellipse.level=0.95, pointsize=4,
			 mean.point=FALSE,
             repel = TRUE # Avoid text overlapping
             )

#clearly, we need to do normalization. center 
# we need to do box plot for each sample
pdf(file="normalization.pdf")      
boxplot(mat)                       
dev.off()

pdf(file="pca_batch.pdf")      
pind                       
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

#now start doing the sample annotation
fullname<-data.frame(FullRunName=names(df)[-c(1:3)])
groups<-cbind(fullname, groups)
groups$MS_batch<-rep(c(1:3),9)
groups$sample_num<-sub(x=groups$samples, pattern="[0-9A-Za-z]+\\-", 
	replacement="")

groups$sample_num<-sub(x=groups$sample_num, pattern="[a-zA-Z]+", 
	replacement="")

groups$sample_num<-sub(x=groups$sample_num, pattern="\\-", 
	replacement="")
groups$sample_num<-as.numeric(groups$sample_num)
groups$group_num<-rep(c(1:3),rep(9,3))

groups$order<-(groups$MS_batch-1)*9+((groups$group_num-1)*3+groups$sample_num)

groups$MS_batch<-as.factor(groups$MS_batch)
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
plot_boxplot(example_proteome, groups,
batch_col = batch_col#, color_scheme = color_list[[batch_col]]
)

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
)

pca2 = plot_PCA(quantile_normalized_matrix, groups, color_by = 'group',
plot_title = 'group'#, color_scheme = color_list[['MS_batch']]
)

library(ggpubr)
ggarrange(pca1, pca2, ncol = 2, nrow = 1)