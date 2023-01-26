rm(list = ls())
library(tidyr)
library(tibble)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(SeuratObject)
logFCfilter=0.5
adjPvalFilter=0.05
levels(sce$RNA_snn_res.0.3)
####################FindMarkers##############
sce.markers <- FindAllMarkers(object = sce,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sce.pos.markers <- FindAllMarkers(object = sce,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = logFCfilter)
sig.markers=sce.markers[(abs(as.numeric(as.vector(sce.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sce.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="markers.xls",sep="\t",row.names=F,quote=F)
sig.pos.markers=sce.pos.markers[(abs(as.numeric(as.vector(sce.pos.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sce.pos.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.pos.markers,file="pos-markers.xls",sep="\t",row.names=F,quote=F)
sce.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top50 <- sce.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.table(top50,file="Top50Genes.txt",sep="\t",col.names= NA)

DefaultAssay(sce) <-"RNA" #"integrated"
Idents(sce)<-sce$seurat_clusters

cluster0.markers <- FindMarkers(sce, ident.1 = 0, min.pct = 0.25,only.pos = TRUE)
head(cluster0.markers, n = 10)

cluster1.markers <- FindMarkers(sce, ident.1 = 1, min.pct = 0.25,only.pos = TRUE)
head(cluster1.markers, n = 10)

cluster2.markers <- FindMarkers(sce, ident.1 = 2, min.pct = 0.25,only.pos = TRUE)
head(cluster2.markers, n = 10)

cluster3.markers <- FindMarkers(sce, ident.1 = 3, min.pct = 0.25,only.pos = TRUE)
head(cluster3.markers, n = 10)

cluster4.markers <- FindMarkers(sce, ident.1 = 4, min.pct = 0.25,only.pos = TRUE)
head(cluster4.markers, n = 10)

cluster5.markers <- FindMarkers(sce, ident.1 = 5, min.pct = 0.25,only.pos = TRUE)
head(cluster5.markers, n = 10)

cluster6.markers <- FindMarkers(sce, ident.1 = 6, min.pct = 0.25,only.pos = TRUE)
head(cluster6.markers, n = 10)

cluster7.markers <- FindMarkers(sce, ident.1 = 7, min.pct = 0.25,only.pos = TRUE)
head(cluster7.markers, n = 10)

cluster8.markers <- FindMarkers(sce, ident.1 = 8, min.pct = 0.25,only.pos = TRUE)
head(cluster8.markers, n = 10)

cluster9.markers <- FindMarkers(sce, ident.1 = 9, min.pct = 0.25,only.pos = TRUE)
head(cluster9.markers, n = 10)

cluster10.markers <- FindMarkers(sce, ident.1 = 10, min.pct = 0.25,only.pos = TRUE)
head(cluster10.markers, n = 10)

cluster11.markers <- FindMarkers(sce, ident.1 = 11, min.pct = 0.25,only.pos = TRUE)
head(cluster11.markers, n = 10)

cluster12.markers <- FindMarkers(sce, ident.1 = 12, min.pct = 0.25,only.pos = TRUE)
head(cluster12.markers, n = 10)

cluster13.markers <- FindMarkers(sce, ident.1 = 13, min.pct = 0.25,only.pos = TRUE)
head(cluster13.markers, n = 10)

cluster14.markers <- FindMarkers(sce, ident.1 = 14, min.pct = 0.25,only.pos = TRUE)
head(cluster14.markers, n = 10)

cluster15.markers <- FindMarkers(sce, ident.1 = 15, min.pct = 0.25,only.pos = TRUE)
head(cluster15.markers, n = 10)

cluster16.markers <- FindMarkers(sce, ident.1 = 15, min.pct = 0.25)
head(cluster16.markers, n = 10)

cluster17.markers <- FindMarkers(sce, ident.1 = 17, min.pct = 0.25)
head(cluster17.markers, n = 10)

cluster18.markers <- FindMarkers(sce, ident.1 = 18, min.pct = 0.25)
head(cluster18.markers, n = 10)

cluster19.markers <- FindMarkers(sce, ident.1 = 19, min.pct = 0.25)
head(cluster19.markers, n = 10)
sce$origCluster<-sce$origCluster[,drop=T]
a<-table(sce$origCluster,sce$ann)%>%as.data.frame()

# cluster2.markers <- FindMarkers(sce, ident.1 = 1,only.pos=TRUE, min.pct = 0.25)
# head(cluster2.markers, n = 10)

# cluster5.markers_1 <- FindMarkers(sce, ident.1 = "Myeloid-CCL5",only.pos=TRUE, min.pct = 0.25)
# head(cluster5.markers_1, n = 10)
# 
# cluster8.markers_1 <- FindMarkers(sce, ident.1 = "Macr-CXCL10",only.pos=TRUE, min.pct = 0.25)
# head(cluster8.markers_1, n = 10)

save(list=ls(),file = "Markers.Rdata")
####################FindConservedMarkers####################################
BiocManager::install('multtest')
install.packages('metap')
library(multtest)
library(metap)
load("D:/Bioinfrolf/data/tmpRdata/01.ALLimmune.combined.Rdata")

head(immune.combined@meta.data)
immune.combined<-sce
DefaultAssay(immune.combined) <- "RNA"
nk.markers0 <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers1 <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers2 <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers3 <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers4 <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers5 <- FindConservedMarkers(immune.combined, ident.1 = 5, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers6 <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers7 <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers8 <- FindConservedMarkers(immune.combined, ident.1 = 8, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers9 <- FindConservedMarkers(immune.combined, ident.1 = 9, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers10 <- FindConservedMarkers(immune.combined, ident.1 = 10, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers11 <- FindConservedMarkers(immune.combined, ident.1 = 11, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers12 <- FindConservedMarkers(immune.combined, ident.1 = 12, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
nk.markers14 <- FindConservedMarkers(immune.combined, ident.1 = 14, grouping.var = "GSE", verbose = FALSE,only.pos = TRUE)
write.csv(nk.markers1,file = "immune_makers1.csv")
write.csv(nk.markers2,file = "immune_makers2.csv")
write.csv(nk.markers3,file = "immune_makers3.csv")
write.csv(nk.markers4,file = "immune_makers4.csv")
write.csv(nk.markers0,file = "immune_makers0.csv")
write.csv(nk.markers5,file = "immune_makers5.csv")
write.csv(nk.markers6,file = "immune_makers6.csv")
write.csv(nk.markers7,file = "immune_makers7.csv")
write.csv(nk.markers8,file = "immune_makers8.csv")
write.csv(nk.markers9,file = "immune_makers9.csv")
write.csv(nk.markers10,file = "immune_makers10.csv")
write.csv(nk.markers11,file = "immune_makers11.csv")
write.csv(nk.markers12,file = "immune_makers12.csv")
write.csv(nk.markers13,file = "immune_makers13.csv")
sce.markers<-rbind(nk.markers0,nk.markers1,nk.markers2,nk.markers3,nk.markers4,
                   nk.markers5,nk.markers6,nk.markers7,nk.markers8,nk.markers9,
                   nk.markers10,nk.markers11,nk.markers12,nk.markers13)
head(sce.markers)
sce.pos.markers<-sce.markers[sce.markers$GSE140228_avg_log2FC>0&sce.markers$GSE156625_avg_log2FC>0,]
write.table(sce.markers,file="markers.xls",sep="\t",row.names=F,quote=F)
write.table(sce.pos.markers,file="pos-markers.xls",sep="\t",row.names=F,quote=F)

clusters<-sce@meta.data$ann%>%as.factor()%>%levels()
sce@meta.data%>%head() 
levels(sce)

############annotation###################
Idents(sce)%>%unique()
immune.combined<-sce
immune.combined$ann%>%unique()
Idents(immune.combined)<-immune.combined$seurat_clusters
immune.combined <- RenameIdents(immune.combined, `0` = "CD8 T cells", `1` = "CD4 T cells", `2` = "NK cells",
                                `3` = "Myeloid cells",`4` = "Tregs", `5` = "HSPs T/NK cells",`6` = "HSPs CD4 T cells",`7` = "Myeloid cells",
                                `8` = "Plasma cells",`9` = "B cells", `10` = "T-NK-Cycle", `11` = "NK cells", 
                                `12` = "Mast cells", `13` = "ILCs")
immune.combined@meta.data$global<-Idents(immune.combined)
Idents(immune.combined)<-immune.combined$seurat_clusters

immune.combined <- RenameIdents(immune.combined, `0` = "T-CD8-GZMK", `1` = "T-CD4-IL7R", `2` = "NK-GNLY",
                                `3` = "Myeloid-S100A9",`4` = "Treg-CTLA4", `5` = "T/NK-HSPA6",`6` = "T-CD4-HSPA6",`7` = "Myeloid-APOC1",
                                `8` = "B-Plasma-MZB1",`9` = "B-CD79A", `10` = "T/NK-Cycle-STMN1", `11` = "NK-TM4SF1", 
                                `12` = "Mast-TPSB2", `13` = "ILC-LILRA4")

immune.combined@meta.data$ann<-Idents(immune.combined)
Idents(sce)<-sce$seurat_clusters
save(list = ls(),file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/01.ALLimmune.combined.Rdata")

Idents(sce)<-sce$seurat_clusters
levels(sce)
sce <- RenameIdents(sce, `0` = "Macr-APOE", `1` = "Macr-RPs", `2` = "DC-CD1C",
                    `3` = "Mono-VCAN",`4` = "Myeloid-CCL5", `5` = "Mono-CD16",`6` = "DC-CD1C-RPs",
                    `7` = "Macr-STAT1",`8` = "Macr-SPP1",`9` = "DC-CLEC9A", `10` = "DC-MKI67",
                    `11` = "DC-LAMP3", `12` = "Macr-VCAM1", `13` = "pDC-IGKC")

save(list = ls(),file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/04.harmony-myeloid.Rdata")

noRPs <- RenameIdents(noRPs, `0` = "Macr-APOE", `1` = "Macr-TMSB4X", `2` = "DC-CD1C",
                    `3` = "Mono-VCAN",`6` = "Myeloid-CCL5", `4` = "Mono-CD16",
                    `5` = "Macr-STAT1",`7` = "Macr-SPP1",`8` = "DC-CLEC9A", `9` = "DC-MKI67",
                    `11` = "DC-LAMP3", `12` = "Macr-VCAM1", `10` = "pDC-IGKC")
noRPs$ann1<-Idents(noRPs)
Idents(noRPs)<-noRPs$seurat_clusters

sce1$ann<-sce$ann
t<-table(noRPs$ann,noRPs$ann1)%>%as.data.frame()
levels(t$Var1)
colnames(t)<-c("ann","noribo","count")
t<-dcast(t,ann~noribo)
rownames(t)<-t$ann
t<-t[,-1]
t1<-t
for(i in 1:nrow(t)){
  for(j in 1:ncol(t)){
    t1[i,j]<-t[i,j]/apply(t, 1, sum)[i]
  }
}
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(0,0.5,1), c("light grey","white", "#FF0000"))
jpeg(file="No-RPs Myeliod Clusters.jpeg",width=5,height=4,units = "in", res = 1000)
Heatmap(
  t1,
  col = col_fun,
  name = 'percentage',
  column_title = 'No-RPs Myeloid Clusters',
  row_title = 'Myeloid Clusters',
  show_row_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12)
)
dev.off()
save(list = ls(),file = "noribo_new.Rdata")
############singleR#############################
library(SingleR)
library(celldex)
library(scater)
load("D:/Bioinfrolf/Bioinfrolf/data/BlueprintEncodeData.Rdata")
head(sce@meta.data$seurat_clusters)
test.seu<-sce
test.count=test.seu@assays$RNA@data
clusters<-test.seu@meta.data$seurat_clusters
common_bp <- intersect(rownames(test.count), rownames(bp.se))
bp.se <- bp.se[common_bp,]
test.count_forbp <- test.count[common_bp,]
test.count_forbp.se <- SummarizedExperiment(assays=list(counts=test.count_forbp))
test.count_forbp.se <- logNormCounts(test.count_forbp.se)
pred.main.bp<-SingleR(test=test.count,ref =bp.se,
                      labels = bp.se$label.main,
                      clusters = clusters,
                      assay.type.test = "logcounts",
                      assay.type.ref = "logcounts")
str(pred.main.bp,max.level = 3)
celltype=data.frame(clusterID=rownames(pred.main.bp),
                    celltype=pred.main.bp$labels,
                    stringsAsFactors = F)

all.markers <- metadata(pred.main.bp)$de.genes
levels(test.seu)
new.cluster.ids <- celltype$celltype
current.cluster.ids<-levels(test.seu)
test.seu@meta.data$SingleRCluster <- plyr::mapvalues(x = test.seu@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

head(test.seu@meta.data)
p1<-plotScoreHeatmap(pred.main.bp)
jpeg(file = "HN-ScoreHeatmap.jpg", width = 5, height = 6, units = "in", res =1000)
plotScoreHeatmap(pred.main.bp,clusters=sce$ann1)
                 #annotation_colors	=list("Clusters"=col14))#,as.data.frame(colData(sce)[,"donor",drop=FALSE]))
dev.off()
test.count_forbp.se$labels <- test.seu@meta.data$SingleRCluster 
test.count_forbp.se$ann<-test.seu@meta.data$ann
test.count_forbp.se$cluster<-Idents(test.seu)

p2<-plotHeatmap(test.count_forbp.se, order_columns_by="ann",column_annotation_colors=list("ann"=col14),
                features=unique(unlist(all.markers$DC))[1:50])
jpeg(file = "ScoreHeatmap-DCs.jpg", width = 5, height = 6, units = "in", res =1000)
p2
dev.off()

p3<-plotHeatmap(test.count_forbp.se, order_columns_by="ann",column_annotation_colors=list("ann"=col14),
                features=unique(unlist(all.markers$Macrophages))[1:50])
jpeg(file = "ScoreHeatmap-Macrophages.jpg", width = 5, height = 6, units = "in", res =1000)
p3
dev.off()

p4<-plotHeatmap(test.count_forbp.se, order_columns_by="ann",column_annotation_colors=list("ann"=col14),
                features=unique(unlist(all.markers$Monocytes))[1:50])
jpeg(file = "ScoreHeatmap-Monocytes.jpg", width = 5, height = 6, units = "in", res =1000)
p4
dev.off()

jpeg(file = "GSE115-UMAP-singleR.jpg", width = 6, height = 4, units = "in", res =1000)
UMAPPlot(test.seu,group.by="SingleRCluster")
dev.off()

save(test.count_forbp.se,pred.main.bp,file = "SingleRforMyeloids.Rdata")
sce<-test.seu
############seleck MNPs####################
sce<- subset(x = sce, 
             subset = ann%in%c("Macr-APOE","Macr-RPs","Mono-VCAN","Mono-CD16",
                                            "Macr-STAT1","Macr-SPP1","Macr-VCAM1"))
sce$ann<-sce$ann[,drop=T]
unique(sce$ann)
current.cluster.ids <-levels(as.factor(sce$tissue_sub))
new.cluster.ids <-c("AN","CT","PT")
sce$tissue_sub<- plyr::mapvalues(x = sce$tissue_sub,
                                 from = current.cluster.ids, to = new.cluster.ids)
save(sce,col8,file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/03.MNPs.Rdata")

