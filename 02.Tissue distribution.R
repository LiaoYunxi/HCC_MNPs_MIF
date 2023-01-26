rm(list = ls())
setwd("D:/Bioinfrolf/HCC-SC/myeloid/figs")
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr_noCCL5.Rdata")
library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(Seurat)
library(ggplot2)
library(ggrepel)

sce<-Mono_Macr
head(sce@meta.data)
table(df$group)
pcadim=23
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

df <- sce@meta.data[,c("tissue","tissue_sub","ann","GSE")]
head(df)
df$ann<-as.character(df$ann)
df$cluster <- unlist(lapply(strsplit(as.character(df$ann),split="-"),function(x){x[1]}))
df$tissue_sub<-factor(df$tissue_sub,levels = c("AN", "PT","CT" ))
df$group<-factor(df$group,levels = c("Healthy","Normal", "TumorEdge","TumorCore" ))
df$ann%>%table()
df$ann<-factor(df$ann,levels = c("Mono-VCAN","Mono-CD16","Macr-VCAM1","Macr-STAT1",
                                 "Macr-APOE","Macr-RPs","Macr-SPP1"))

#df <- df[!df$cluster%in%c('DC'),]

# df <- df[grep("CD14CD16",df$cluster_cancer,invert=T),]
# df <- df[df$tissue!='L',]

jpeg(file="Proportion-MNNPs-TissueSub.jpeg",width=5,height=4,units = "in", res = 1000)
ggplot(df, aes(factor(ann)))+ 
  geom_bar(aes(fill = tissue_sub), position = "fill")+ xlab("")+
  ylab("Proportion")+  theme_bw() +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
  )+scale_y_continuous(labels = scales::label_percent())+
  #labs(title ="Liver MNPs",x = '',y = "percent")+
  scale_fill_manual(values=c("#00A08A","#F98400","#FF0000" ))
dev.off()

summary <- table(df[,c('ann','tissue_sub')])
roe <- as.data.frame(ROIE(summary))

summary <- table(df[,c('ann','tissue')])
roe <- as.data.frame(ROIE(summary))
#roe$marker <- rownames(roe)
#roe$cancer <- unlist(lapply(strsplit(roe$marker,"-"),function(x){x[3]}))

class(roe)
roe1<-as.matrix(roe)
col <- colorRampPalette(c("red", "white", "blue"))(256)
# roe1<-data.frame(Nomal=roe$Normal,Tumor=roe$Tumor)
# rownames(roe1)<-rownames(roe)
jpeg(file="ROE-Tissue.jpeg",width=6,height=5,units = "in", res = 1000)
pheatmap::pheatmap(roe,scale = "none",
                   color = colorRampPalette(c("#FEFEC0","#D73027"))(100),
                   #annotation_col = group,
                   cellwidth = 40, cellheight = 30,
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = FALSE, cluster_col = FALSE,
                   #display_numbers=T,
                   border_color ='white',
                   number_format="%.2f",number_color="black",
                   angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                   main ="Tissue distribution (Ro/e)")
dev.off()


pheatmap::pheatmap(roe,scale = "none",
                   color = colorRampPalette(c("#FEFEC0","#D73027"))(100),
                   #annotation_col = group,
                   cellwidth = 40, cellheight = 40,
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = FALSE, cluster_col = FALSE,
                   display_numbers=T,
                   border_color ='white',
                   number_format="%.2f",number_color="black",
                   angle_col=0,fontsize=15,fontsize_col=12,fontsize_row=12,
                   main ="Tissue distribution (Ro/e)")
roe<-read.table("roe.txt")
