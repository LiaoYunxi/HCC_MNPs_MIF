rm(list = ls())
options(stringsAsFactors = F)
library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(stringr)
library(harmony)(Seurat)
library(ggplot2)
library(cowplot)
library(SeuratObject)
library(clustree)
#############zhang########################
setwd("D:/Bioinfrolf/Bioinfrolf/data")
Droplet_cellinfo<-read.table("GSE140228_UMI_counts_Droplet_cellinfo.tsv",sep="\t",header=T,check.names=F)
str(Droplet_cellinfo)
dim(Droplet_cellinfo)
sum(as.character(Droplet_cellinfo$celltype_global))
Droplet_cellinfo[1:5,]
Droplet_barcodes<-read.table("GSE140228_UMI_counts_Droplet_barcodes.tsv",sep="\t",header=F,check.names=F)
str(Droplet_barcodes)
head(Droplet_barcodes)
dim(Droplet_barcodes)
Droplet_genes<-read.table("GSE140228_UMI_counts_Droplet_genes.tsv",sep="\t",header=T,check.names=F)
str(Droplet_genes)
Droplet_genes[1:5,1:5]
dim(Droplet_genes)
Droplet <- Matrix::readMM("GSE140228_UMI_counts_Droplet.mtx")
str(Droplet)
class(Droplet)#dgTMatrixϡ??????
dim(Droplet)
Droplet[1:5,1:5]
levels(as.factor(Droplet_cellinfo$celltype_global))
levels(as.factor(Droplet_cellinfo$celltype_sub))
levels(as.factor(Droplet_cellinfo$Tissue))
levels(as.factor(Droplet_cellinfo$Tissue_sub))

table(as.factor(Droplet_cellinfo$celltype_global),as.factor(Droplet_cellinfo$celltype_sub))
cell_all<-Droplet_cellinfo[Droplet_cellinfo$Tissue%in%c("Normal","Tumor"),]
cell_table<-table(as.factor(cell_all$celltype_global),as.factor(cell_all$celltype_sub))%>%as.data.frame()
cell_table<-cell_table[cell_table$Freq!="0",]
colnames(Droplet)<-Droplet_cellinfo$Barcode
rownames(Droplet)<-Droplet_genes$SYMBOL
im_matrix<-Droplet[,colnames(Droplet)%in%cell_all$Barcode]
dim(im_matrix)
save(im_matrix,cell_all,Droplet_genes,file = "D:/Bioinfrolf/immune sc/ztmp.Rdata")
rm(Droplet,Droplet_barcodes,Droplet_cellinfo)

#############sanma##########################
setwd("D:/Bioinfrolf/data")
setwd("D:/Bioinfrolf/single cell/GSE156625/supply/GSE156625F")
HCC_barcodes<-read.table("GSE156625_HCCFbarcodes.tsv",sep="\t",header=F,check.names=F)
dim(HCC_barcodes)
HCC_featutes<-read.table("GSE156625_HCCFgenes.tsv",sep="\t",header=F,check.names=F)
dim(HCC_featutes)
HCC_matrix<-Matrix::readMM("GSE156625_HCCFmatrix.mtx")
dim(HCC_matrix)
str(HCC_matrix)
HCC_matrix[1:5,1:5]
head(HCC_barcodes)
head(HCC_featutes)
colnames(HCC_matrix)<-HCC_barcodes$V1
rownames(HCC_matrix)<-HCC_featutes$V2

library(plyr)
library(dplyr)
library(scales)
library(reticulate)
options(stringsAsFactors=FALSE)

parse_h5ad <- function(adata){
  require(reticulate)
  ad <- import("anndata", convert = FALSE)
  ada <- ad$read_h5ad(adata)
  meta <- py_to_r(ada$obs)
  if(class(ada$raw$X)[1] == "scipy.sparse.csr.csr_matrix"){
    exp <- t(py_to_r(ada$raw$X$toarray()))
  }
  else{
    exp <- t(py_to_r(ada$raw$X))
  }
  rownames(exp) <- rownames(py_to_r(ada$raw$var))
  colnames(exp) <- rownames(meta)
  return(
    list(
      metadata = meta,
      expression = exp
    )
  )
}

h5ad <- parse_h5ad("GSE156625_HCCFscanpyobj.h5ad")

# py_config()
Sys.which("python")
#activate  C:\Users\ADMINI~1\AppData\Local\R-MINI~1\envs\r-reticulate
#use_python("C:/ProgramData/Anaconda3/python.exe")
#use_virtualenv("base(root)")

head(h5ad$metadata)
cluster<-read.table("cluster-all.txt",sep="\t",header=F,check.names=F)
head(cluster)
index<-c()
for (i in 1:nrow(cluster)){
  a<-unlist(strsplit(cluster$V1[i]," ",fixed = TRUE))[1]
  index<-c(index,a)
}

str(h5ad)
h5ad$metadata<-h5ad$metadata[h5ad$metadata$louvain%in%index,]
h5ad$expression[1:5,1:5]
h5ad$expression<-h5ad$expression[,rownames(h5ad$metadata)]

index<-c()
for (i in 1:nrow(cluster)){
  a<-paste(unlist(strsplit(cluster$V1[i]," ",fixed = TRUE))[-1],sep = " ",collapse =" ")
  index<-c(index,a)
}

h5ad$metadata$louvain<-h5ad$metadata$louvain[,drop = TRUE]
new.cluster.ids <-index
current.cluster.ids <-levels(h5ad$metadata$louvain)

#h5ad$metadata$Cluster <- sapply(as.vector(h5ad$metadata$MajorCluster), function(x){unlist(strsplit(x,"_"))[2]})
h5ad$metadata$Cluster <- plyr::mapvalues(x = h5ad$metadata$louvain, from = current.cluster.ids, to = new.cluster.ids)
head(h5ad$metadata)
levels(as.factor(h5ad$metadata$patient_tumorsection))
levels(as.factor(h5ad$metadata$PIC))
levels(as.factor(h5ad$metadata$PNC))
table(as.factor(h5ad$metadata$PNC),as.factor(h5ad$metadata$PIC))
table(as.factor(h5ad$metadata$NormalvsTumor),as.factor(h5ad$metadata$PNC))
table(h5ad$metadata$Cluster,as.factor(h5ad$metadata$PIC))
head(h5ad$metadata$PNC)
save(h5ad,file = "h5ad.Rdata")
metas<-h5ad$metadata
save(metas,file = "mata.GSE156625.Rdata")
{
  HCC_matrix[1:5,1:5]
  meta<-h5ad$metadata
  HCC_matrix<-HCC_matrix[,rownames(h5ad$metadata[h5ad$metadata$Cluster%in%index,])]
  dim(HCC_matrix)
  meta1<-meta[colnames(HCC_matrix),]
  core_matrix<-HCC_matrix[,rownames(meta1[meta1$PNC=="C",])]
  para_matrix<-HCC_matrix[,rownames(meta1[meta1$PNC=="P",])]
  normal_matrix<-HCC_matrix[,rownames(meta1[meta1$PNC=="N",])]
  
  save(core_matrix,para_matrix,normal_matrix,meta,file = "PNC-noTNK.Rdata")
}

#############DoubletFinder#####################
setwd("D:/Bioinfrolf/Bioinfrolf/tmpRdata")
load("00.GSE156625h5ad.Rdata")
im_matrix1<-h5ad$expression
rm(h5ad)
load("00.GSE140228tmp.Rdata")
sce <- CreateSeuratObject(counts = im_matrix1, min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = sce@var.genes, npcs = 50, verbose = FALSE)

sce[["percent.mt"]] <- PercentageFeatureSet(object = sce, pattern = "^MT-")
sce<- subset(x = sce, subset = nFeature_RNA > 200 & percent.mt < 10)

ElbowPlot(sce,ndims=50,reduction = "pca")
pcadim<-30

sce <- sce %>% 
  RunUMAP(dims = 1:pcadim) %>% 
  FindNeighbors( dims = 1:pcadim) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

library(DoubletFinder)
sweep.res.list <- paramSweep_v3(sce, PCs = 1:pcadim, sct = FALSE,num.cores=1)
for(i in 1:length(sweep.res.list)){
  if(length(sweep.res.list[[i]]$pANN[is.nan(sweep.res.list[[i]]$pANN)]) != 0){
    if(i != 1){
      sweep.res.list[[i]] <- sweep.res.list[[i - 1]]
    }else{
      sweep.res.list[[i]] <- sweep.res.list[[i + 1]]
    }
  }
}

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
nExp_poi <- round(0.03*length(colnames(sce)))
sce <- doubletFinder_v3(sce, 
                        PCs = 1:pcadim, pN = 0.25, pK = pk_good, 
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

colnames(sce@meta.data)[ncol(sce@meta.data)]="DoubletFinder"
head(sce@meta.data)
levels(as.factor(sce$DoubletFinder))
table(sce$DoubletFinder)


jpeg(file = "umap-Doublets.jpg", width = 6, height = 4, units = "in", res = 2000)
DimPlot(sce,reduction = "umap",pt.size = 0.03,
        group.by = "DoubletFinder",cols =wes_palettes$FantasticFox1[c(5,3)])
dev.off()

DimPlot(sce,reduction = "umap",pt.size = 0.03,group.by = "DoubletFinder0.05")
write.csv(sce@meta.data[,7:8],file = 'doublets_s.csv')

sce<-sce[,rownames(sce@meta.data[sce$DoubletFinder=="Singlet",])]

#############CCA##################
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/00.ALLmatrix.Rdata")
DoubleListS<-read.csv('D:/Bioinfrolf/Bioinfrolf/table/doublets_s.csv')
DoubleListZ<-read.csv('D:/Bioinfrolf/Bioinfrolf/table/doublets_z.csv')
table(DoubleListZ$DoubletFinder)
DoubleListS<-DoubleListS[DoubleListS$DoubletFinder=="Doublet",]
DoubleListZ<-DoubleListZ[DoubleListZ$DoubletFinder=="Doublet",]
doublets<-c(DoubleListZ$X,DoubleListS$X)
im_matrix<-im_matrix[,!colnames(im_matrix)%in%c(DoubleListZ$X)]
dim(im_matrix)
im_matrix1<-im_matrix1[,!colnames(im_matrix)%in%c(DoubleListS$X)]

load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/01.ALLmatrix.Rdata")
zm<-CreateSeuratObject(counts = im_matrix, min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000)

sm<-CreateSeuratObject(counts = im_matrix1, min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000)

features <- SelectIntegrationFeatures(object.list = c(zm,sm))
immune.anchors <- FindIntegrationAnchors(object.list = c(zm,sm), anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
rm(im_matrix,im_matrix1)

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)

immune.combined@meta.data$GSE<-ifelse(rownames(immune.combined@meta.data)%in%colnames(zm),"GSE140228","GSE156625")
table(immune.combined$GSE)
#####################Plot#############################
setwd("D:/HCC-SC/all immune/figs")

jpeg(file = "PCA.jpg", width = 25, height = 30, units = "in", res =300)
VizDimLoadings(object = immune.combined, dims = 1:30,nfeatures = 20,ncol = 5)
dev.off()

jpeg(file = "PCAheatmap.jpg", width = 20, height = 30, units = "in", res =300)
DimHeatmap(object = immune.combined, dims = 1:30, cells = 500, balanced = TRUE,nfeatures = 20,ncol=5)
dev.off()

# immune.combined <- JackStraw(object = immune.combined, num.replicate = 100,dims=50)
# immune.combined <- ScoreJackStraw(object = immune.combined, dims = 1:50)
# pdf(file="pcaJackStraw.pdf",width=8,height=6)
# JackStrawPlot(object = immune.combined, dims = 1:30)
# dev.off()

jpeg(file = "ElbowPlot.jpg", width = 5, height = 4, units = "in", res =1000)
ElbowPlot(immune.combined,ndims=50)
dev.off()
pcadim<-30

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:pcadim)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:pcadim)
#immune.combined <- FindClusters(immune.combined, resolution = 0.5)
immune.combined <- FindClusters(immune.combined, resolution = c(seq(.1,1.6,.1)))
head(immune.combined@meta.data)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:pcadim)

jpeg(file = "clustree.jpg", width = 10, height = 12, units = "in", res =1000)
clustree(immune.combined@meta.data, prefix = "integrated_snn_res.")
dev.off()

p1=DimPlot(immune.combined, reduction = 'umap', group.by = 'integrated_snn_res.0.5',
           label = TRUE, pt.size = 0.05) + NoLegend()
p2=DimPlot(immune.combined, reduction = 'umap',group.by = 'integrated_snn_res.0.4',
           label = TRUE, pt.size = 0.05) + NoLegend()
library(patchwork)
pdf(file="resolution-02vs03.pdf",width=8,height=4)
p1+p2
dev.off()

head(immune.combined@meta.data)
head(immune.anchors@object.list)
zv<-immune.anchors@object.list[[1]]%>%VariableFeatures()
sv<-immune.anchors@object.list[[2]]%>%VariableFeatures()
match(zv,sv)%>%na.omit()
dim(immune.combined@meta.data)
immune.combined@meta.data$GSE<- c(rep("GSE140228", 32113), rep("GSE156625", 44540))

#immune.combined<-FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 4000)
#change the Idents
immune.combined@meta.data$seurat_clusters<-immune.combined@meta.data$integrated_snn_res.0.2
Idents(immune.combined)
levels(immune.combined)
head(immune.combined@meta.data)
new.cluster.ids <- immune.combined@meta.data$integrated_snn_res.0.2
names(new.cluster.ids) <- rownames(immune.combined@meta.data)
Idents(immune.combined)<-new.cluster.ids
str(immune.combined@assays$integrated)
immune.combined[["percent.mt"]] <- PercentageFeatureSet(object = immune.combined,pattern = "^MT-",assay = "RNA")
immune.combined[["percent.ribo"]] <- PercentageFeatureSet(object = immune.combined,pattern = "^RP[SL]",assay = "RNA")
jpeg(file = "featureViolin.jpg", width = 24, height =5, units = "in", res =300)
VlnPlot(object = immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
dev.off()

immune.combined<- subset(x = immune.combined, subset = nFeature_RNA > 200 & percent.mt < 10)

jpeg(file="featureCor.jpg",width=12,height=8,units = "in", res =1000) 
plot1 <- FeatureScatter(object = immune.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        pt.size=0.2,group.by = "GSE",cols = col2)
plot2 <- FeatureScatter(object = immune.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        pt.size=0.2,group.by = "GSE",cols = col2)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# top10 <- head(x = VariableFeatures(object = immune.combined), 10)
# pdf(file="featureVar.pdf",width=10,height=6)  
# plot1 <- VariableFeaturePlot(object = immune.combined)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
# dev.off()
table(immune.combined$ann,immune.combined$origAnn)
############celltype#################
origAnn<-c(cell_all$celltype_global,cell_all1$Cluster)
origCluster<-c(cell_all$celltype_sub,cell_all1$louvain)
ann<-c(cell_all$celltype_sub,cell_all1$Cluster)
barcodes<-c(cell_all$Barcode,rownames(cell_all1))
tissue<-c(cell_all$Tissue,cell_all1$NormalvsTumor)
tissue_sub<-c(cell_all$Tissue_sub,cell_all1$PNC)
donor<-c(cell_all$Donor,cell_all1$patientno)
table(cell_all$Sample)
table(cell_all$Tissue_sub)
metadata<-data.frame(origAnn=origAnn,origCluster=origCluster,ann=ann,tissue=tissue,donor=donor,tissue_sub=tissue_sub)
rownames(metadata)<-barcodes
metadata<-metadata[rownames(immune.combined@meta.data),]

immune.combined@meta.data$tissue_sub<-metadata$tissue_sub
immune.combined@meta.data$tissue<-metadata$tissue
immune.combined@meta.data$origAnn<-metadata$origAnn
immune.combined@meta.data$origCluster<-metadata$origCluster
immune.combined@meta.data$ann<-metadata$ann
table(immune.combined$tissue_sub)
current.cluster.ids <-levels(as.factor(immune.combined@meta.data$tissue_sub))
current.cluster.ids
new.cluster.ids <-c("TumorCore","Normal","TumorEdge","Normal","TumorCore","TumorEdge")
immune.combined@meta.data$tissue_sub <- plyr::mapvalues(x = immune.combined@meta.data$tissue_sub,
                                                        from = current.cluster.ids, to = new.cluster.ids)

current.cluster.ids <-levels(as.factor(immune.combined@meta.data$tissue))
current.cluster.ids
new.cluster.ids <-c("Normal","Tumor","Normal","Tumor")
immune.combined@meta.data$tissue <- plyr::mapvalues(x = immune.combined@meta.data$tissue,
                                                    from = current.cluster.ids, to = new.cluster.ids)
levels(as.factor(immune.combined@meta.data$ann))
levels(as.factor(cell_all1$Cluster))
levels(as.factor(immune.combined@meta.data$origCluster))
current.cluster.ids <-levels(as.factor(immune.combined@meta.data$ann))
current.cluster.ids
new.cluster.ids <-c(levels(as.factor(cell_all1$Cluster)),current.cluster.ids[8:length(current.cluster.ids)])
immune.combined@meta.data$ann <- plyr::mapvalues(x = immune.combined@meta.data$ann,
                                                 from = current.cluster.ids, to = new.cluster.ids)
celltype_table<-as.matrix(table(as.factor(immune.combined@meta.data$ann),as.factor(immune.combined@meta.data$seurat_clusters)))
celltype_table1<-as.matrix(table(as.factor(immune.combined@meta.data$origAnn),as.factor(immune.combined@meta.data$seurat_clusters)))
write.csv (celltype_table1, file ="celltype1.csv")
write.csv (celltype_table, file ="celltype.csv")

pdf(file="umap-TvsN.pdf",width=10,height=4)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = .1, split.by = "tissue")
dev.off()
save(immune.combined,file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/01.immune.combined.Rdata")


#immune.combined <- RunTSNE(immune.combined , dims = 1:pcadim)
# pdf(file="TSNE.pdf",width=6,height=4)
# TSNEPlot(object = immune.combined,reduction = "tsne",group.by = "seurat_clusters",pt.size = .1, label = TRUE)  
# dev.off()
