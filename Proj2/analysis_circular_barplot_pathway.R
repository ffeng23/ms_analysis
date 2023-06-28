#R code to do circular barplot
# for all pathways
# ref to do string wrap!!!
#https://stackoverflow.com/questions/21878974/wrap-long-axis-labels-via-labeller-label-wrap-in-ggplot2
#

#
library(here)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(ggpubr)

#read the pathways
#xt.down<-read_excel(
#	here("Proj2","VISTA_CD4+T pathways analysis_innateDB_05262023.xlsx"), 
#	sheet="InnateDB_Pathway_VISTA-Down"
#	)
xt.down<-read_excel(
  here("Proj2","Signaling pathway_disease groups.xlsx"), 
  sheet=1
  )
xt.down.clean <- xt.down %>% 
	dplyr::filter(`Pathway p-value (corrected)`< 0.05)
#xt.down.clean<-xt.down[xt.down$'Pathway p-value (corrected)'<0.05,]
xt.down.clean <- xt.down.clean %>%
  mutate(group= replace(Catalog, is.na(Catalog),"Others")) 


#now get up and down
upDown<-which(duplicated(xt.down.clean$`Pathway Name`))

xt.down.clean.upDown<-xt.down.clean[upDown,]
xt.down.clean.upDown<-as.data.frame(xt.down.clean.upDown)

xt.down.clean<-xt.down.clean[-upDown,]

xt.down.clean<-as.data.frame(xt.down.clean)
rownames(xt.down.clean)<-xt.down.clean$`Pathway Name`
#set the up and down
xt.down.clean[xt.down.clean.upDown$`Pathway Name`,10]<-"Up&Down"

#455 rows

xt.down.clean<-  xt.down.clean %>%
	dplyr::rename(p.adj=`Pathway p-value (corrected)`,
		p=`Pathway p-value`, Genes=`Genes (Symbol|IDBG-ID|Ensembl|Entrez|Fold Change|P-Value)`,
		pathway.count=`Genes in InnateDB for this entity`) %>%
	select(`Pathway Name`, `Pathway Id`,p.adj, pathway.count, Catalog,Genes,group) %>%
	mutate(p.adj=-log10(p.adj)) 

xt.down.clean$group<-paste0(xt.down.clean$group,"(", xt.down.clean$Genes,")")
#xt.down.clean$group<-paste0("(", xt.down.clean$Genes,")",xt.down.clean$group)



#make up some group for now <------------#
#xt.down.clean$group<-1
#xt.down.clean$group[1:(floor(dim(xt.down.clean)[1]/5)*5)]<-
#	rep(seq(1:5),floor(dim(xt.down.clean)[1]/5))

xt.down.clean$group<-as.factor(xt.down.clean$group)

data<-xt.down.clean
#data<-data[order(data$Genes),]
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
to_add$Genes<-to_add$group
to_add$Genes <- sub(x=to_add$Gene,pattern=".+\\(([a-zA-Z&]+)\\)",replacement="\\1")
data <- rbind(data, to_add)
data <- data %>% arrange(Genes,group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label, 
# label individual pathway names
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for group names.
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))



# prepare a data frame for grid (scales), at the end of each group
grid_data <- base_data
#grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
#grid_data$start <- grid_data$start - 1
#grid_data <- grid_data[-1,]
grid_data$start<-grid_data$end+empty_bar
grid_data$end<-grid_data$end+1

max.value<-max(xt.down.clean$p.adj)
min.value<-min(xt.down.clean$p.adj)
data$value<-data$p.adj/max.value *100
data$id<-factor(data$id, levels=data$id)
p2 <- ggplot(data, aes(x=id, y=value, fill=Genes)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=id, y=value, fill=Genes), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  #geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend =40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("","","","") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  #geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-80,105) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-4,4), "cm"),
    plot.title = element_text(hjust = 0.3, vjust = -60,face="bold") 
  ) +
  coord_polar() + #ggtitle("Down-/Up-regulated Pathways")+
  #geom_text(data=label_data, aes(x=id, y=p.adj+10, label='Pathway Name', hjust=hjust), 
  #		color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -15, label=group),
  	hjust=0.5,#c(1,1,0,0,0), 
    colour = "black", alpha=0.8, size=2, 
  	fontface="bold", inherit.aes = FALSE)

p <- ggplot(data, aes(x=id, y=value, fill=Genes)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=id, y=value, fill=Genes), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  #geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend =40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("","","","") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  #geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-20,105) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-4,4), "cm"),
    plot.title = element_text(hjust = 0.3, vjust = -60,face="bold") 
  ) +
  coord_polar() + #ggtitle("Down-/Up-regulated Pathways")+
  #geom_text(data=label_data, aes(x=id, y=p.adj+10, label='Pathway Name', hjust=hjust), 
  #   color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  
  


  pdf(file=here("Proj2/removeReps/"
    ,"vista_updown_circular.pdf"), width=7,height=7)
  p2
  dev.off()

  pdf(file=here("Proj2/removeReps/"
    ,"vista_updown_circular_noname.pdf"), width=7,height=7)
  p
  dev.off()


###section using the old style with separated
#  might not work!!! be careful.

  ###########doing the upregulated 
#read the pathways
xt.up<-read_excel(
	here("Proj2","VISTA_CD4+T pathways analysis_innateDB_05262023.xlsx"), 
	sheet="InnateDB_Pathway_VISTA-Up"
	)
xt.up.clean <- xt.up %>% 
	dplyr::filter(`Pathway p-value (corrected)`< 0.05)
#xt.down.clean<-xt.down[xt.down$'Pathway p-value (corrected)'<0.05,]

#455 rows

xt.up.clean<-  xt.up.clean %>%
	dplyr::rename(p.adj=`Pathway p-value (corrected)`,
		p=`Pathway p-value`,
		pathway.count=`Genes in InnateDB for this entity`) %>%
	select(`Pathway Name`, `Pathway Id`,p.adj, pathway.count) %>%
	mutate(p.adj=-log10(p.adj)) 

#make up some group for now <------------#
xt.up.clean$group<-1
xt.up.clean$group[1:(floor(dim(xt.up.clean)[1]/5)*5)]<-
	rep(seq(1:5),floor(dim(xt.up.clean)[1]/5))

xt.up.clean$group<-as.factor(xt.up.clean$group)

data<-xt.up.clean
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
#grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
#grid_data$start <- grid_data$start - 1
#grid_data <- grid_data[-1,]
grid_data$start<-grid_data$end+empty_bar
grid_data$end<-grid_data$end+1

max.value<-max(xt.up.clean$p.adj)
min.value<-min(xt.up.clean$p.adj)
data$value<-data$p.adj/max.value *100
p2 <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend =40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("","","","") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  #geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-80,105) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-4,4), "cm"),
    plot.title = element_text(hjust = 0.3, vjust = -60,face="bold") 
  ) +
  #scale_y_continuous(minor_breaks = c(-80, 0,2,4,6))+
  coord_polar() + ggtitle("Up-regulated Pathways")+
  geom_text(data=label_data, aes(x=id, y=p.adj/max.value*100+2, label=`Pathway Name`, hjust=hjust), 
  		color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -15, label=group),
  	hjust=c(1,1,0,0,0), colour = "black", alpha=0.8, size=4, 
  	fontface="bold", inherit.aes = FALSE)

  pdf(file=here(output.dir,"vista_up_circular.pdf"), width=7,height=7)
  p2
  dev.off()

