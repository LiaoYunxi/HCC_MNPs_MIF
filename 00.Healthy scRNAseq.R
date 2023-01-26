options(stringsAsFactors = F)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(stringr)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(ggrepel)
library(SeuratObject)
library(Seurat)
library(harmony)
library(clustree)
library(edgeR)
library(limma)
setwd("D:/Bioinfrolf/new figs/")
rm(list=ls())
####################GSE115469####
GSE115469_Data<-read.csv("D:/Bioinfrolf/Bioinfrolf/data/GSE115469_Data.csv",header=T,check.names=F)
rownames(GSE115469_Data)<-GSE115469_Data[,1]
GSE115469_Data<-GSE115469_Data[,-1]
GSE115469_Data[1:5,1:5]
GSE115469matrix<-apply(GSE115469_Data, 1, function(x){
  expm1(x)*10000
})

GSE115<-CreateSeuratObject(count=GSE115469_Data)
dim(GSE115@assays$RNA@counts)
GSE115@assays$RNA@data<-GSE115@assays$RNA@counts
GSE115<-FindVariableFeatures(GSE115,selection.method = "vst", nfeatures = 3000)
GSE115<-ScaleData(GSE115,verbose = FALSE) %>% 
  RunPCA(pc.genes = GSE115@var.genes, npcs = 50, verbose = FALSE)

GSE115[["percent.mt"]] <- PercentageFeatureSet(object = GSE115, pattern = "^MT-")

jpeg(file="GSE115-featureViolin.jpg",width=15,height=4,units = "in", res = 1000)
VlnPlot(object = GSE115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
GSE115$sample<-GSE115$orig.ident

GSE115<- subset(x = GSE115, subset = percent.mt < 10)
#GSE115<- subset(x = GSE115, subset = nFeature_RNA > 200)

jpeg(file="GSE115mye-featureCor.jpg",width=10,height=4,units = "in", res = 1000) 
plot1 <- FeatureScatter(object = GSE115, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5,group.by = "sample")
plot2 <- FeatureScatter(object = GSE115, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5,group.by = "sample")
CombinePlots(plots = list(plot1, plot2))
dev.off()

top10 <- head(x = VariableFeatures(object = GSE115), 10)

plot <- VariableFeaturePlot(object = GSE115)
label.data <- plot$data[top10, ]
label.data$labels <- top10
plot1 <- plot + geom_text_repel(aes(x = mean,y = variance.standardized,label = labels), data = label.data)
#CombinePlots(plots = list(plot1, plot2))
jpeg(file="GSE115mye-featureVar.jpg",width=6,height=4,units = "in", res = 1000)
plot1
dev.off()

jpeg(file="GSE115mye-pcaGene.jpeg",width = 25, height = 30, units = "in", res =300)
VizDimLoadings(object = GSE115, dims = 1:30,nfeatures = 20,ncol = 5,reduction = "pca")
dev.off()

jpeg(file="GSE115mye-PCA.jpeg",width=6,height=4,units = "in", res = 1000)
DimPlot(object = GSE115, reduction = "pca",group.by = "sample")
dev.off()

jpeg(file="GSE115mye-pcaHeatmap.jpeg",width=20,height=30,units = "in", res = 300)
DimHeatmap(object = GSE115, dims = 1:30, cells = 500, balanced = TRUE,nfeatures = 20,ncol=5,reduction = "pca")
dev.off()

GSE115 <- JackStraw(object = GSE115, num.replicate = 100,dims=50)
GSE115 <- ScoreJackStraw(object = GSE115, dims = 1:50)
jpeg(file="pcaJackStraw.jpeg",width=6,height=4,units = "in", res = 1000)
JackStrawPlot(object = GSE115, dims = 1:30)
dev.off()

jpeg(file="GSE115mye-ElbowPlot.jpeg",width = 5, height = 4, units = "in", res =1000)
ElbowPlot(GSE115,ndims=50,reduction = "harmony")
dev.off()
pcadim<-22

GSE115 <- GSE115 %>% 
  RunUMAP(reduction = "pca", dims = 1:pcadim) %>% 
  FindNeighbors(reduction = "pca", dims = 1:pcadim) %>% 
  FindClusters(resolution = c(seq(.3,0.8,.1))) %>% 
  identity()

jpeg(file="GSE115mye-clustree.jpg",width=10,height=12,units = "in", res = 1000)
clustree(GSE115@meta.data, prefix = "RNA_snn_res.")
dev.off()
colnames(GSE115@meta.data)

p1=DimPlot(GSE115, reduction = 'umap', group.by = 'RNA_snn_res.0.3',
           label = TRUE, pt.size = 0.1) + NoLegend()
p2=DimPlot(GSE115, reduction = 'umap',group.by = 'SingleRCluster',
           label = TRUE, pt.size = 0.1) + NoLegend()
library(patchwork)
jpeg(file="GSE115-resolution-03vs05.jpeg",width=8,height=4,units = "in", res = 1000)
p1+p2
dev.off()

jpeg(file="GSE115-feature-TNK.jpeg",width =10,height = 8,units = "in", res = 1000)
FeaturePlot(GSE115, features = c(
  "CD8A","IL7R","FOXP3","GNLY"),reduction = 'umap',
  cols = wes_palettes$FantasticFox1[c(3,5)], min.cutoff = "q4")
dev.off()

jpeg(file="GSE115-feature-myeloid.jpeg",width =10,height = 12,units = "in", res = 1000)
FeaturePlot(GSE115, features = c(
  "LYZ","CD163","HLA-DRA","CD68","VCAM1","MARCO"),reduction = 'umap',
  cols = wes_palettes$FantasticFox1[c(3,5)], min.cutoff = "q6")
dev.off()

sce$seurat_clusters<-sce$RNA_snn_res.0.3
sce <- RenameIdents(sce, `0` ="0", `1` = "Monocytes", `2` = "CD8 T cells",`3` = "3",
                    `4` = "Plasma cells",`5` = "NK cells",`6` = "NK cells",`7` = "B cells",`8` = "Hepatocytes",
                    `9` = "Red blood cell",`10` = "cycle T cells")
sce <- RenameIdents(sce, `0` = "Macr-RPs",`1` = "Macr-APOE", `2` = "Mono-VCAN",
                    `3` = "Mono-VCAN", `4` = "Kupffer cells")
sce$ann<-Idents(sce)
Idents(sce)<-sce$seurat_clusters
save(list = ls(),file = "GSE115469sce.Rdata")

myeloidGSE115<-rownames(GSE115@meta.data[GSE115$seurat_clusters%in%c(1),])
GSE115<-GSE115[,myeloidGSE115]

GSE115mye<-sce
save(GSE115mye,file = "GSE115469mye.Rdata")

GSE115469_Data<-GSE115469_Data[,monoGSE115]
write.table(GSE115469_Data,file = "myeloidGSE115469.txt",quote = F)
GSE115<-GSE115[,myeloidGSE115]

mnGSE115<-rownames(GSE115@meta.data[GSE115$RNA_snn_res.0.3%in%c(0,1,4),])
GSE115469_Data<-GSE115469_Data[,mnGSE115]
write.table(GSE115469_Data,file = "mnGSE115469.txt",quote = F)
GSE115<-GSE115[,mnGSE115]
save(GSE115,file = "GSE115469MN.Rdata")

GSE115mye$seurat_clusters<-GSE115mye$RNA_snn_res.0.3
Idents(GSE115mye)<-GSE115mye$seurat_clusters

####################GSE129933####
GSE129933metadata<-read.table("D:/Bioinfrolf/single cell/fudan/GSE129933_cell_metadata.tsv",sep="\t",header=T,check.names=F)
GSE129933matrix<-read.table("D:/Bioinfrolf/single cell/fudan/GSE129933_count_matrix.tsv",sep="\t",header=T,check.names=F)
rownames(GSE129933matrix)<-GSE129933matrix$gene
GSE129933matrix<-GSE129933matrix[,-1]
table(GSE129933metadata$sample,GSE129933metadata$group)
GSE129933metadata<-GSE129933metadata[GSE129933metadata$group=="Non-Diseased",]
GSE129933matrix<-GSE129933matrix[,GSE129933metadata$cell]
rownames(GSE129933metadata)<-GSE129933metadata$cell
dim(GSE129933matrix)
GSE129933matrix[1:5,1:5]
table(GSE129933metadata$cluster)

GSE129<-CreateSeuratObject(count=GSE129933matrix,meta.data = GSE129933metadata[,-1])%>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)
head(GSE129@meta.data)
GSE129[["percent.mt"]] <- PercentageFeatureSet(object = GSE129, pattern = "^MT-")
GSE129[["percent.ribo"]] <- PercentageFeatureSet(object = GSE129, pattern = "^RP[SL]")

jpeg(file="GSE129-featureViolin.jpg",width=15,height=4,units = "in", res = 1000)
VlnPlot(object = GSE129, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="GSE129-featureCor.jpg",width=10,height=4,units = "in", res = 1000) 
plot1 <- FeatureScatter(object = GSE129, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5,group.by = "sample")
plot2 <- FeatureScatter(object = GSE129, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5,group.by = "sample")
CombinePlots(plots = list(plot1, plot2))
dev.off()

top10 <- head(x = VariableFeatures(object = GSE129), 10)

plot <- VariableFeaturePlot(object = GSE129)
label.data <- plot$data[top10, ]
label.data$labels <- top10
plot1 <- plot + geom_text_repel(aes(x = mean,y = variance.standardized,label = labels), data = label.data)
jpeg(file="GSE129-featureVar.jpg",width=6,height=4,units = "in", res = 1000)
plot1
dev.off()

jpeg(file="GSE129-pcaGene.jpeg",width = 25, height = 30, units = "in", res =300)
VizDimLoadings(object = GSE129, dims = 1:30,nfeatures = 20,ncol = 5,reduction = "pca")
dev.off()

jpeg(file="GSE129-PCA.jpeg",width=6,height=4,units = "in", res = 1000)
DimPlot(object = GSE129, reduction = "pca",group.by = "sample")
dev.off()

jpeg(file="GSE129-pcaHeatmap.jpeg",width=20,height=30,units = "in", res = 300)
DimHeatmap(object = GSE129, dims = 1:30, cells = 500, balanced = TRUE,nfeatures = 20,ncol=5,reduction = "pca")
dev.off()

GSE129 <- JackStraw(object = GSE129, num.replicate = 100,dims=50)
GSE129 <- ScoreJackStraw(object = GSE129, dims = 1:50)
jpeg(file="pcaJackStraw.jpeg",width=6,height=4,units = "in", res = 1000)
JackStrawPlot(object = GSE129, dims = 1:30)
dev.off()

jpeg(file="GSE129-ElbowPlot.jpeg",width = 5, height = 4, units = "in", res =1000)
ElbowPlot(GSE129,ndims=50,reduction = "pca")
dev.off()
pcadim<-30

GSE129 <- GSE129 %>% 
  RunUMAP(reduction = "pca", dims = 1:pcadim) %>% 
  FindNeighbors(reduction = "pca", dims = 1:pcadim) %>% 
  FindClusters(resolution = c(seq(.3,0.8,.1))) %>% 
  identity()

jpeg(file="GSE129-clustree.jpg",width=10,height=12,units = "in", res = 1000)
clustree(GSE129@meta.data, prefix = "RNA_snn_res.")
dev.off()
colnames(GSE129@meta.data)

p1=DimPlot(GSE129, reduction = 'umap', group.by = 'RNA_snn_res.0.8',
           label = TRUE, pt.size = 0.1) + NoLegend()
p2=DimPlot(GSE129, reduction = 'umap',group.by = 'SingleRCluster',
           label = TRUE, pt.size = 0.1) + NoLegend()
library(patchwork)
jpeg(file="GSE129-resolution-03vs05.jpeg",width=8,height=4,units = "in", res = 1000)
p1+p2
dev.off()

jpeg(file="GSE129-feature-TNK.jpeg",width =10,height = 8,units = "in", res = 1000)
FeaturePlot(GSE129, features = c(
  "CD8A","IL7R","FOXP3","GNLY"),reduction = 'umap',
  cols = wes_palettes$FantasticFox1[c(3,5)], min.cutoff = "q4")
dev.off()

jpeg(file="GSE129-feature-myeloid.jpeg",width =10,height = 12,units = "in", res = 1000)
FeaturePlot(GSE129, features = c(
  "LYZ","CD163","HLA-DRA","CD68","VCAM1","MARCO"),reduction = 'umap',
  cols = wes_palettes$FantasticFox1[c(3,5)], min.cutoff = "q6")
dev.off()

save(GSE129,file = "tmp.Rdata")
sce$seurat_clusters<-sce$RNA_snn_res.0.3
table(sce$cluster,sce$seurat_clusters)
sce <- RenameIdents(sce, `0` ="Lymphocytes", `1` = "Fibroblasts&LEC&PEC", `2` = "Monocytes",
                    `3` = "B cells",`4` = "Erythrocytes")
sce$ann<-Idents(sce)
Idents(sce)<-sce$seurat_clusters
save(list = ls(),file = "GSE129933sce.Rdata")
p1<-UMAPPlot(sce,group.by="ann",pt.size=3)
p2<-UMAPPlot(sce,group.by="cluster",pt.size=3)
p1+p2

GSE129mye<-GSE129[,rownames(GSE129@meta.data[GSE129$ann=="Monocytes",])]
save(GSE129mye,file = "GSE129933mye.Rdata")
# GSE129933metadata<-GSE129933metadata[GSE129933metadata$cluster=="Monocytes",]
# rownames(GSE129933matrix)<-GSE129933matrix$gene
# GSE129933matrix<-GSE129933matrix[,-1]
# GSE129933matrix[1:5,1:5]
# GSE129933matrix<-GSE129933matrix[,GSE129933metadata$cell]
# 
# GSE129933_Data<-NormalizeData(GSE129933matrix)%>%as.matrix()#sce@assays$RNA@data%>%as.matrix()
# GSE129933_Data[1:5,1:5]
# dim(GSE129933_Data)
# rm(GSE129933matrix)
####################GSE136103############
setwd("D:/Bioinfrolf/Bioinfrolf/data")
if(F){
  Droplet_barcodes<-read.table("GSM4041150_healthy1_cd45+_barcodes.tsv",sep="\t",header=F,check.names=F)
  Droplet_genes<-read.table("GSM4041150_healthy1_cd45+_genes.tsv",sep="\t",header=F,check.names=F)
  Droplet <- Matrix::readMM("GSM4041150_healthy1_cd45+_matrix.mtx")
  dim(Droplet)
  nrow(Droplet)
  Droplet_barcodes$V1<-unlist(lapply(Droplet_barcodes$V1,function(x){
    paste0("hn1_",x)
  }))
  rownames(Droplet)<-Droplet_genes$V2
  colnames(Droplet)<-Droplet_barcodes$V1
  d1<-Droplet
  
  Droplet_barcodes<-read.table("GSM4041153_healthy2_cd45+_barcodes.tsv",sep="\t",header=F,check.names=F)
  Droplet_genes<-read.table("GSM4041153_healthy2_cd45+_genes.tsv",sep="\t",header=F,check.names=F)
  Droplet <- Matrix::readMM("GSM4041153_healthy2_cd45+_matrix.mtx")
  Droplet_barcodes$V1<-unlist(lapply(Droplet_barcodes$V1,function(x){
    paste0("hn2_",x)
  }))
  rownames(Droplet)<-Droplet_genes$V2
  colnames(Droplet)<-Droplet_barcodes$V1
  d2<-Droplet
  
  Droplet_barcodes<-read.table("GSM4041155_healthy3_cd45+_barcodes.tsv",sep="\t",header=F,check.names=F)
  Droplet_genes<-read.table("GSM4041155_healthy3_cd45+_genes.tsv",sep="\t",header=F,check.names=F)
  Droplet <- Matrix::readMM("GSM4041155_healthy3_cd45+_matrix.mtx")
  Droplet_barcodes$V1<-unlist(lapply(Droplet_barcodes$V1,function(x){
    paste0("hn3_",x)
  }))
  rownames(Droplet)<-Droplet_genes$V2
  colnames(Droplet)<-Droplet_barcodes$V1
  d3<-Droplet
  
  Droplet_barcodes<-read.table("GSM4041158_healthy4_cd45+_barcodes.tsv",sep="\t",header=F,check.names=F)
  Droplet_genes<-read.table("GSM4041158_healthy4_cd45+_genes.tsv",sep="\t",header=F,check.names=F)
  Droplet <- Matrix::readMM("GSM4041158_healthy4_cd45+_matrix.mtx")
  Droplet_barcodes$V1<-unlist(lapply(Droplet_barcodes$V1,function(x){
    paste0("hn4_",x)
  }))
  rownames(Droplet)<-Droplet_genes$V2
  colnames(Droplet)<-Droplet_barcodes$V1
  d4<-Droplet
  rm(Droplet,Droplet_barcodes,Droplet_genes)
}

sce <- CreateSeuratObject(counts = d1, min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)
sce1<-sce
h1<-rownames(sce1@meta.data[sce1$DoubletFinder=="Singlet",])

index<-c(h1,h2,h3,h4)#
sce <- CreateSeuratObject(counts = cbind(d1,d2,d3,d4)[,index], min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)

save(d1,d2,d3,d4,h1,h2,h3,h4,sce1,sce2,sce3,sce4,file = "GSE136103matrix.Rdata")
sce@meta.data$sample<- c(rep("h1",length(h1)), rep("h2", length(h2)), rep("h3",length(h3)), rep("h4", length(h4)))

{
  setwd("D:/Bioinfrolf/new figs/GSE115/")  
  options(repr.plot.height = 4, repr.plot.width = 10)
  p1 <- DimPlot(object = sce, reduction = "pca", pt.size = .1, group.by = "sample")
  p2 <- VlnPlot(object = sce, features = "PC_1", group.by = "sample", pt.size = .1)
  jpeg(file="GSE115-PCA.jpg",width=10,height=4,units = "in", res = 1000)
  plot_grid(p1,p2)
  dev.off()
}

#options(repr.plot.height = 2.5, repr.plot.width = 6)

sce <- sce %>% 
  RunHarmony("sample", plot_convergence = TRUE,max.iter.harmony = 50,max.iter.cluster = 20)
harmony_embeddings <- Embeddings(sce, 'harmony')
harmony_embeddings[1:5, 1:5]
pcadim<-23

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = sce, reduction = "harmony",group.by = "sample")
p2 <- VlnPlot(object = sce, features = "harmony_1",group.by = "sample")
jpeg(file="harmony_GSE136.jpg",width=10,height=4,units = "in", res = 1000)
plot_grid(p1,p2)
dev.off()

p1=DimPlot(sce, reduction = 'umap', group.by = 'RNA_snn_res.0.3',
           label = TRUE, pt.size = 0.1) + NoLegend()
p2=DimPlot(sce, reduction = 'umap',group.by = 'SingleRCluster',
           label = TRUE, pt.size = 0.1) + NoLegend()
library(patchwork)
jpeg(file="GSE115-SingleRCluster.jpg",width=10,height=4,units = "in", res = 1000)
p1+p2
dev.off()

sce <- RenameIdents(sce, `0` ="CD8 T cells", `1` = "CD8 T cells", `2` = "CD4 T cells",`3` = "CD4 T cells",
                    `4` = "NK cells",`5` = "Monocytes",`6` = "B cells",`7` = "Kupffer cells",`8` = "T/Nk cells",
                    `9` = "DCs",`10` = "cycle T cells",'11'="Macrophages","12"="Plasma cells",
                    '13'="B cells","14"="Mast cells","15"="Hepatocytes")

sce$seurat_clusters<-sce$RNA_snn_res.0.3
Idents(sce)<-sce$seurat_clusters
sce$ann<-Idents(sce)
Idents(sce)<-sce$seurat_clusters
table(sce$seurat_clusters)
save(list = ls(),file = "GSE136sce.Rdata")
GSE136mye<-subset(x = sce, subset =ann%in%c("Monocytes","Macrophages","Kupffer cells","DCs"))
save(GSE136mye,file = "GSE136mye.Rdata")

rm(mono)
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
load("D:/Bioinfrolf/single cell/fudan/GSE115469-liver/GSE115469sce.Rdata")
monoGSE115<-read.table(file = "monoGSE115.txt")
GSE115<-GSE115[,monoGSE115$x]

meta1<-Mono_Macr@meta.data
data1<-Mono_Macr@assays$RNA@data%>%as.matrix()
dim(data)
meta2<-GSE136@meta.data
data2<-GSE136@assays$RNA@data%>%as.matrix()
meta3<-GSE115@meta.data
data3<-GSE115@assays$RNA@data%>%as.matrix()
data4<-GSE129933_Data
meta4<-GSE129933metadata
gene<-Reduce(intersect,list(rownames(data1),rownames(data2),rownames(data3)))
data1<-data1[gene,]
data2<-data2[gene,]
data3<-data3[gene,]
gene<-Reduce(intersect,list(rownames(data4),rownames(data2),rownames(data3)))
data4<-data4[gene,]
save(data1,data2,data3,data4,meta1,meta2,meta3,meta4,file = "data.Rdata")
save(data2,data3,data4,meta2,meta3,meta4,file = "data-hn.Rdata")
####################HN#################
gene<-intersect(rownames(GSE115@assays$RNA@counts),rownames(GSE136mye@assays$RNA@counts))
gene<-intersect(gene,rownames(GSE129mye@assays$RNA@counts))
im_matrix<-GSE115@assays$RNA@data[gene,]
im_matrix1<-GSE136mye@assays$RNA@data[gene,]
im_matrix2<-GSE129mye@assays$RNA@data[gene,]

sce<-CreateSeuratObject(count=cbind(im_matrix,im_matrix1), min.cells = 5)
sce@assays$RNA@data<-sce@assays$RNA@counts
sce<-FindVariableFeatures(sce,selection.method = "vst", nfeatures = 3000)
sce<-ScaleData(sce,verbose = FALSE) 
sce@assays$RNA@data<-sce@assays$RNA@counts
sce<-RunPCA(sce,pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)
sce@meta.data$GSE<- c(rep("GSE115469", ncol(im_matrix)), rep("GSE136103", ncol(im_matrix1)))
pcadim<-18
sce <- RenameIdents(sce, `0` ="Kupffer cells", `2` = "Mono-VCAN", `3` = "CD4 T cells",`3` = "CD4 T cells",
                    `4` = "NK cells",`5` = "Monocytes",`6` = "Macrophages",`7` = "Kupffer cells",`8` = "T/Nk cells",
                    `9` = "DCs",`10` = "cycle T cells",'11'="Macrophages","12"="Plasma cells",
                    '13'="B cells","14"="Mast cells","15"="Hepatocytes")

features <- SelectIntegrationFeatures(object.list = c(GSE115mye,GSE136mn),
                                      nfeatures = 3000,fvf.nfeatures = 3000,)
immune.anchors <- FindIntegrationAnchors(object.list =c(GSE115mye,GSE136mn), anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors,normalization.method = c("LogNormalize"))
DefaultAssay(immune.combined) <- "integrated"
immune.combined@meta.data$GSE<-ifelse(rownames(immune.combined@meta.data)%in%colnames(GSE115mye),"GSE115469","GSE136103")
setwd("D:/Bioinfrolf/new figs/HN")
pcadim<-20