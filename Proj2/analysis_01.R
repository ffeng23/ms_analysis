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