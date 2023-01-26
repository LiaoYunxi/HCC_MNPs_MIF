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
load("slingshot.Rdata")
sce@colData$ann<-Mono_Macr@meta.data$ann
sce <- slingshot::slingshot(sce, clusterLabels = 'ann', reducedDim = 'PCA')

summary(sce$slingPseudotime_1)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

library(tradeSeq)
# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
#heatdata <- assays(sce)$norm[topgenes, pst.ord]
heatclus <- sce$ann[pst.ord]

class(log1p(heatdata))
jpeg(file="D:/HCC-SC/myeloid/TI method/slingshot-heatmap1-1.jpg",width =10,height = 8,units = "in", res = 2000)
heatmap(as.matrix(log1p(heatdata)), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
dev.off()

jpeg(file="D:/HCC-SC/myeloid/TI method/slingshot-heatmap2-1.jpg",width =10,height = 8,units = "in", res = 2000)
heatmap(as.matrix(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
dev.off()

sce@colData$GMM%>%as.factor()%>%levels()
table(sce@colData$GMM,sce@colData$ann)

if(plot){
  lin1 <- getLineages(rd1, cl1, start.clus = '1')
  lin1
  plot(rd1, col =col8[cl1], asp = 1, pch = 16)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  
  lin2 <- getLineages(rd2, colData(sce)$ann, start.clus = 'Mono-VCAN')
  lin1
  plot(rd2, col = col8[colData (sce) $ ann], asp = 1, pch = 16)
  lines(SlingshotDataSet(lin2), lwd = 3, col = 'black')
  lin2 <- getLineages(rd1, cl1, start.clus= '4', end.clus = '1')
  lin2 <- getLineages(rd2, colData(sce)$ann, start.clus = 'Mono-VCAN',end.clus = 'Macr-SPP1') 
  lin2
  
  
  plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
  lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
  
  plot(rd1, col = brewer.pal(9,"Set1")[sce$ann], asp = 0.5, pch = 16)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  levels(Mono_Macr)
  
  rd1%>%head()
  d1<-as.data.frame(rd1)
  
  col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6","#E2D200","#B40F20","#273046","#FD6467")
  ggplot(d1, aes(x = PC1, y = PC2, color = sce$ann)) + 
    geom_point(size = 0.5, alpha = 1)+
    labs(title ="Myeloid cell",
         x = 'PC1',
         y = "PC1")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 15,face = 'bold'),
          axis.text = element_text(size=12), # Font size of axis labels.
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
          #legend.position = c(.92, .9),                     # Position of the legend.
          #legend.margin = margin(6, 6, 6, 6),                # Margin of the legend.
          #legend.background = element_rect(size = 0.2, colour = 1)
    )+
    scale_color_manual(values=brewer.pal(8,"Set1"))
  
  save(list = ls(),file = "slingshot-2.Rdata")
  write.table(topgenes,file = "slingshot-topgene.txt",row.names = F,col.names = F,quote = F)
}