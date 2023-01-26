rm(list = ls())
options(stringsAsFactors = F)
library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(stringr)
library(harmony)
library(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratObject)
library(clustree)
library(ggrepel)
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
Mono_Macr<-subset(x = Mono_Macr,subset = ann!="Macr-CCL5")
Mono_Macr@meta.data$ann<-Mono_Macr@meta.data$ann[,drop=T]
sce1<-Mono_Macr
rm(Mono_Macr)
###################myeloids##################
setwd("D:/Bioinfrolf/data/tmpRdata")
load("h5ad.Rdata")
load("ztmp.Rdata")
cell_all1<-h5ad$metadata
im_matrix1<-h5ad$expression
rm(Droplet_genes,h5ad)
im_matrix<-Matrix::as.matrix(im_matrix)
head(cell_all1)
gene<-intersect(rownames(im_matrix),rownames(im_matrix1))
length(gene)
im_matrix1<-im_matrix1[gene,]
im_matrix<-im_matrix[gene,]

setwd("D:/Bioinfrolf/data/All")
rm(list = ls())
load("01.matrix.Rdata")
cell_mye<-cell_all[cell_all$celltype_global%in%c("Myeloid","Myeloid-Mast","Myeloid-Liver-doublets"),]
cell_mye<-cell_all[cell_all$celltype_global%in%c("Myeloid"),]
cell_mye$celltype_global<-cell_mye$celltype_global[,drop=T]
levels(as.factor(cell_mye$celltype_sub))
levels(as.factor(cell_all1$Cluster))

cell_mye1<-cell_all1[cell_all1$Cluster%in%c("Myeloid"),]
cell_mye1<-cell_all1[cell_all1$Cluster%in%c("Myeloid","Mast Cells"),]
cell_mye1$Cluster<-cell_mye1$Cluster[,drop=T]
levels(as.factor((cell_mye1$Cluster)))
im_matrix<-im_matrix[,cell_mye$Barcode]
im_matrix1<-im_matrix1[,rownames(cell_mye1)]
save(cell_mye,cell_mye1,im_matrix,im_matrix1,file = "04.Myeloidmatrix.Rdata")
###################harmony##########################
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/02.Myeloidmatrix.Rdata")
sce <- CreateSeuratObject(counts = cbind(im_matrix,im_matrix1), min.cells = 5)%>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)
sce@meta.data$GSE<- c(rep("GSE140228", ncol(im_matrix)), rep("GSE156625", ncol(im_matrix1)))
sce@meta.data$orig.ident<-sce@meta.data$GSE
sce[["percent.mt"]] <- PercentageFeatureSet(object = sce, pattern = "^MT-")
sce[["percent.ribo"]] <- PercentageFeatureSet(object = sce, pattern = "^RP[SL]")
features = NULL
ribo <- features %||% grep(pattern = "^RP[SL]", x = rownames(sce@assays$RNA@counts),value = TRUE)
mt<- features %||% grep(pattern = "^MT-", x = rownames(sce@assays$RNA@counts),value = TRUE)
sce<- subset(x = sce, subset = percent.mt < 10)

counts <- GetAssayData(object = sce, slot = "counts")
filtered_counts <- counts[!rownames(counts)%in%c(ribo,mt),]
dim(counts)
dim(filtered_counts)
# nonzero <- counts > 0
# keep_genes <- Matrix::rowSums(nonzero) >= 0
# filtered_counts <- counts[keep_genes, ]
sce <- CreateSeuratObject(filtered_counts, meta.data = sce@meta.data)%>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)

# countLOW=quantile(sce$nCount_RNA, prob=c(0.005))  
# countHIGH=quantile(sce$nCount_RNA, prob=0.995) 
# sce <- subset(sce, subset = nCount_RNA > countLOW & nCount_RNA < countHIGH )
featureLOW=quantile(sce$nFeature_RNA, prob=0.01)
sce<- subset(x = sce, subset = nFeature_RNA > featureLOW)

head(sce@meta.data)
sce <- sce %>% RunHarmony("GSE", plot_convergence = TRUE,max.iter.harmony = 50,max.iter.cluster = 20)
harmony_embeddings <- Embeddings(sce, 'harmony')
harmony_embeddings[1:5, 1:5]

VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), group.by = "GSE",
        pt.size=NULL,ncol = 4)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA"), group.by = "GSE",
        pt.size=NULL,ncol = 2)
###################plot#######################
# pdf(file="harmony-featureCor.pdf",width=20,height=6) 
# plot1 <- FeatureScatter(object = sce, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
# plot2 <- FeatureScatter(object = sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
# CombinePlots(plots = list(plot1, plot2))
# dev.off()

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = sce, reduction = "harmony", pt.size = .1, group.by = "GSE",cols = col2)
p2 <- VlnPlot(object = sce, features = "harmony_1", group.by = "GSE", pt.size = .1,cols = col2)
jpeg(file="harmony_GSE-myeloid.jpg",width=10,height=4,units = "in", res = 1000)
plot_grid(p1,p2)
dev.off()

top10 <- head(x = VariableFeatures(object = sce), 10)

plot<- VariableFeaturePlot(object = sce)
label.data <- plot$data[top10, ]
label.data$labels <- top10
plot1 <- plot + geom_text_repel(aes(x = mean,y = variance.standardized,label = labels), data = label.data)
#CombinePlots(plots = list(plot1, plot2))
jpeg(file="harmony-featureVar.jpg",width=6,height=4,units = "in", res = 1000)
plot1
dev.off()

jpeg(file="harmony-pcaGene.jpeg",width = 25, height = 30, units = "in", res =300)
VizDimLoadings(object = sce, dims = 1:30,nfeatures = 20,ncol = 5,reduction = "harmony")
dev.off()

jpeg(file = "harmony-PCAheatmap.jpg", width = 20, height = 30, units = "in", res =300)
DimHeatmap(object = sce, dims = 1:30, cells = 500, balanced = TRUE,nfeatures = 20,ncol=5,reduction = "harmony")
dev.off()

# sce <- JackStraw(object = sce, num.replicate = 100,dims=50)
# sce <- ScoreJackStraw(object = sce, dims = 1:50)
# pdf(file="pcaJackStraw.pdf",width=8,height=6)
# JackStrawPlot(object = sce, dims = 1:30)
# dev.off()

jpeg(file = "harmony-ElbowPlot.jpg", width = 5, height = 4, units = "in", res =1000)
ElbowPlot(sce,ndims=50,reduction = "harmony")
dev.off()
pcadim<-25
###################clustering####################################
set.seed(1)
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = 1:pcadim) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pcadim) %>% 
  FindClusters(resolution = c(seq(.3,0.6,.1))) %>% 
  identity()

jpeg(file = "harmony-clustree.jpg", width = 10, height = 12, units = "in", res =1000)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
dev.off()

colnames(sce@meta.data)

p1=DimPlot(sce, reduction = 'umap', group.by = 'RNA_snn_res.0.3',
           label = TRUE, pt.size = 0.1) + NoLegend()
p2=DimPlot(sce, reduction = 'umap',group.by = 'RNA_snn_res.0.4',
           label = TRUE, pt.size = 0.1) + NoLegend()
library(patchwork)
pdf(file="harmony-resolution-03vs05.pdf",width=8,height=4)
p1+p2
dev.off()

#change the Idents
sce@meta.data$seurat_clusters<-sce@meta.data$RNA_snn_res.0.4
Idents(sce)<-sce$seurat_clusters

jpeg(file="umap_GSE-myeloid.jpg",width=10,height=4,units = "in", res = 1000)
DimPlot(sce, reduction = "umap", group.by = "GSE", pt.size = .1, split.by = 'GSE')
dev.off()

jpeg(file="umap_GSE-myeloid1.jpg",width=6,height=4,units = "in", res = 1000)
DimPlot(sce, reduction = "tsne", group.by = "GSE", pt.size = .1)
dev.off()

jpeg(file="umap.jpeg",width=6,height=4,units = "in", res = 1000)
DimPlot(sce1, reduction = "umap",group.by = "seurat_clusters",label = TRUE, pt.size = .1)
dev.off()

#sce <- RunTSNE(sce , dims = 1:pcadim)
pdf(file="TSNE.pdf",width=6,height=4)
TSNEPlot(object = sce,reduction = "tsne",group.by = "seurat_clusters",pt.size = .1, label = TRUE)  
dev.off()

############meta celltype#################
origCluster<-c(cell_mye$celltype_sub,cell_mye1$louvain)
origAnn<-c(cell_mye$celltype_sub,cell_mye1$Cluster)
barcodes<-c(cell_mye$Barcode,rownames(cell_mye1))
tissue<-c(cell_mye$Tissue,cell_mye1$NormalvsTumor)
donor<-c(cell_mye$Donor,cell_mye1$patientno)
metadata<-data.frame(origAnn=origAnn,origCluster=origCluster,tissue=tissue,donor=donor)
rownames(metadata)<-barcodes
metadata<-metadata[rownames(sce@meta.data),]

sce$donor<-metadata[colnames(sce),]$donor

h<-sce[,colnames(sce)%in%rownames(cell_all1[cell_all1$patientno=="HN1",])]
dh<-h@assays$RNA@data[rownames(sce),]


sce@meta.data$donor<-metadata$donor
sce@meta.data$tissue<-metadata$tissue
sce@meta.data$origAnn<-metadata$origAnn
sce@meta.data$origCluster<-metadata$origCluster
current.cluster.ids <-levels(as.factor(sce@meta.data$tissue))
current.cluster.ids
new.cluster.ids <-c("Normal","Tumor","Normal","Tumor")
sce@meta.data$tissue <- plyr::mapvalues(x = sce@meta.data$tissue,
                                        from = current.cluster.ids, to = new.cluster.ids)

pdf(file="umap-TvsN.pdf",width=10,height=4)
DimPlot(sce, reduction = "umap", label = TRUE, pt.size = .1, split.by = "tissue")
dev.off()
save(metadata,sce,file = "allSeurat.Rdata")