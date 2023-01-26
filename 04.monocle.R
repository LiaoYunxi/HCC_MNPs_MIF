if(preprocessing){
  rm(list = ls())
  setwd("D:/Bioinfrolf/data/HAR-MYE/")
  load("Mono_Macr-tmp.Rdata")
  index<-rownames(Mono_Macr@meta.data[Mono_Macr@meta.data$ann %in%c("Macr-C1QA", "Macr-APOC1", "Mono-S100A8","Mono-FCGR3A","Macr-CXCL10","Macr-SPP1","Macr-VCAM1"),])
  Mono_Macr<-Mono_Macr[,index]
  Myeloid<-Mono_Macr
  levels(Myeloid)
  setwd("D:/Bioinfrolf/data/HAR-MYE/monocle")
  monocle.matrix=as.matrix(Myeloid@assays$RNA@data)
  monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
  write.table(monocle.matrix,file="07.monocleMatrix.txt",quote=F,sep="\t",row.names=F)
  monocle.sample=as.matrix(Myeloid@meta.data)
  monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
  head(monocle.sample)
  write.table(monocle.sample,file="07.monocleSample.txt",quote=F,sep="\t",row.names=F)
  
  monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
  monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
  write.table(monocle.geneAnn,file="07.monocleGene.txt",quote=F,sep="\t",row.names=F)
  
  levels(Myeloid)
  levels(as.factor(Myeloid@meta.data$seurat_clusters))%>%na.omit()
  logFCfilter=0.5
  adjPvalFilter=0.05
  library(Seurat)
  Myeloid.markers <- FindAllMarkers(object = Myeloid,
                                    only.pos = FALSE,
                                    min.pct = 0.25,
                                    logfc.threshold = logFCfilter)
  sig.markers=Myeloid.markers[(abs(as.numeric(as.vector(Myeloid.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(Myeloid.markers$p_val_adj))<adjPvalFilter),]
  write.table(sig.markers,file="07.monocleMarkers.txt",sep="\t",row.names=F,quote=F)
}
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

load("cds.Rdata")
if(F){
  cds <- setOrderingFilter(cds, marker$gene)
  dim(cds)
  diff_test_res <- differentialGeneTest(cds,
                                        fullModelFormulaStr = "~Cluster")
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  save(ordering_genes,file = 'ordering_genes_by_Cluster.Rdata')
  load(file = 'ordering_genes_by_Cluster.Rdata')
  head(ordering_genes)
  cds <- setOrderingFilter(cds, ordering_genes)
}

cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)

library(monocle)
cds <- reduceDimension(cds, max_components = 3,verbose = F,reduction_method = 'tSNE')
cds <- orderCells(cds)
save(cds,file = "cds-DDR-marker.Rdata")
load("cds-DDR-marker.Rdata")

detailed_cell_type_color <- c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
                              "Macr-CCL5"="#B40F20","Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
                              "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A")
col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6","#E2D200","#B40F20","#273046","#FD6467")

cds$State

options(repr.plot.width=3, repr.plot.height=4)
plot_complex_cell_trajectory(cds, color_by = 'State', show_branch_points = T,
                             cell_size = 2, cell_link_size = 0.3, root_states = c(1))

jpeg(file="DDRtree-marker-tree.jpeg",width =12,height = 8,units = "in", res = 2000)
plot_complex_cell_trajectory(cds, color_by = 'cell_type2', show_branch_points = T,
                             cell_size = 1, cell_link_size = 0.5, root_states = c(1)) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme (legend.position="left", legend.title=element_blank()) + scale_color_manual(values = detailed_cell_type_color)
dev.off()


pdf(file="cellType.trajectory2.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
#cell_size = 0.5,cell_link_size = 0.25)
dev.off()

plot_multiple_branches_heatmap(cds[c(positive_score_genes, negtive_score_genes),],
                               branches=c(1, 3, 4, 6, 11, 9),
                               branches_name=c("Mono-VCAN", "Mono-FCGR3A", "Macr-CXCL10", "Macr-CCL5", "Macr-APOE", "Macr-SPP1"),
                               show_rownames=T,
                               num_clusters=4)

#############monocle3###############
rm(list = ls())
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
setwd("D:/HCC-SC/myeloid/TI method")
load("D:/Bioinfrolf/data/tmpRdata/05.Mono_Macr-tmp.Rdata")
scRNAsub<-Mono_Macr

set.seed(123) #random number generator
ind <- sample(2, ncol(Mono_Macr), replace = TRUE, prob = c(1000/ncol(Mono_Macr), (1-1000/ncol(Mono_Macr))))
scRNAsub <- Mono_Macr[,ind==1]

data <- GetAssayData(scRNAsub, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNAsub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="ann") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNAsub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="ann") + ggtitle('int.umap')

?cluster_cells()
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)

cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p

# cds <- order_cells(cds) 存在bug，使用辅助线选择root细胞
p + geom_vline(xintercept = seq(4,5,0.25)) + geom_hline(yintercept = seq(-2,-1,0.25))
embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
embed <- subset(embed, UMAP_1 > 4.5 & UMAP_1 < 5 & UMAP_2 > -2.25 & UMAP_2 < -1.75)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE)

##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="predicted.id", 
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
##寻找共表达模块
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
