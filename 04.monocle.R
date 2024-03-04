#############monocle#####################
rm(list = ls())
library(monocle)
library(tidyverse)
setwd("D:/HCC-SC/myeloid/TI method")
monocle.matrix=read.table("07.monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("07.monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("07.monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("07.monocleMarkers.txt",sep="\t",header=T,check.names=F)

head(monocle.geneAnn)
monocle.sample$seurat_clusters%>%as.factor()%>%levels()
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)

names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

clusterRt=read.table("07.monocleClusterAnn.txt",header=F,sep="\t",check.names=F)
clusterAnn=as.character(clusterRt[,2])
names(clusterAnn)=paste0("cluster",clusterRt[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

cds$cell_type2%>%as.factor%>%levels()

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
save(cds,file = "cds.Rdata")
