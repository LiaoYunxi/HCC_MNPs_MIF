#############slingshot################
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(wesanderson)
library("RColorBrewer")
rm(list = ls())
setwd("D:/HCC-SC/myeloid/TI method/slingshot")
load("D:/Bioinfrolf/data/tmpRdata/05.Mono_Macr-tmp.Rdata")
data("slingshotExample")
slingshotExample$rd%>%head()
slingshotExample$cl%>%head()
class(slingshotExample)

setwd("D:/Bioinfrolf/data/HAR-MYE/slingshot")

table()
test.seu<-Mono_Macr
test.count=test.seu@assays$RNA@counts
test.data=test.seu@assays$RNA@data#GetAssayData(test.seu,slot = "data")#鏍囧噯鍖栧悗鐨?
clusters<-test.seu@meta.data$seurat_clusters
sce<- SingleCellExperiment(assays = List(counts = test.count))
rm(test.seu)
sce<- SingleCellExperiment(assays = List(counts = test.count,norm=test.data))

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 50) >= 10
})
sce <- sce[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

assays(sce)$norm <- FQnorm(assays(sce)$counts)
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
#assays(sce)$norm <- assays(sce)$norm
# pca<-Mono_Macr@reductions$pca
# rd1 <- pca@cell.embeddings[,1:2]
# rd2<-Mono_Macr@reductions$umap@cell.embeddings
# colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')


plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
SingleCellExperiment::reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)
cl1 <- Idents(Mono_Macr)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

?Mclust()
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
cl2 <- kmeans(rd1, centers = 8)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

save(list = ls(),file = "slingshot.Rdata")
