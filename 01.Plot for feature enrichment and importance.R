library("org.Hs.eg.db")
library(enrichplot)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(clusterProfiler)
library(Seurat)
library(ggsci)
library("RColorBrewer")
library(wesanderson)
library(scales)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)
library(ggnewscale)
rm(list = ls())
setwd("D:/Bioinfrolf/Bioinfrolf/GSEA/")
immunomodulator<-read.table("D:/Bioinfrolf/Bioinfrolf/GSEA/immunomodulator.txt",header = F)
colnames(immunomodulator)<-c("protein","type","gene")
Genes<-read.table("D:/Bioinfrolf/Bioinfrolf/model/impMatrix.txt",header = T)
genes=Genes$Feature
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=data.frame(symbol=genes,entrezID=entrezIDs)
gene=out$entrezID

#################feature genes enrichment########################
KEGG<-enrichKEGG(gene = G2,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 minGSSize = 10,
                 maxGSSize = 500,
                 qvalueCutoff = 0.05,
                 use_internal_data = FALSE)
go.BP=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05, 
               ont="BP", readable =T)
head(go.BP@result)

showNum=30
if(nrow(go.BP)<showNum){
  showNum=nrow(go.BP)
}

colorSel="qvalue"
head(go.BP@result)
write.table(go.BP@result,file = "feature_GO.BP.BH.txt",quote = F)
write.table(go.BP@result,file = "KEGG.txt",quote = F)

jpeg(file="feature_GO.BP.jpeg",width =10,height =6,units = "in", res = 2000)
enrichplot::dotplot(go.BP, showCategory = showNum, 
                    color = colorSel,orderBy="GeneRatio")
dev.off()

test<-go.BP@result
library(stringr)
gr1 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,2])
test$fold <- (gr1/gr2)/(bg1/bg2)
test$GeneRatio <- (gr1/gr2)
head(test)

test<-test[order(test$fold,decreasing = T),]  
df<-test[30:1,]
df$Description = factor(df$Description,
                        levels=df$Description)
jpeg(file="feature_GO.BP.fold.jpeg",width =10,height =6,units = "in", res = 2000)
ggplot(df,aes(x = fold,y = Description))+
  geom_point(aes(color = p.adjust,size = Count))+
  scale_color_gradient(low ="#FF0000", high = "#00A08A")+
  xlab("Fold Enrichment")+theme_bw()+
  guides(color = guide_colorbar(reverse = TRUE))
dev.off()

test<-test[order(test$GeneRatio,decreasing = T),] 
df<-test[30:1,]
df$Description = factor(df$Description,
                        levels=df$Description)
jpeg(file="feature_GO.BP.GeneRatio.jpeg",width =10,height =6,units = "in", res = 2000)
ggplot(df,aes(x = GeneRatio,y = Description))+
  geom_point(aes(color = p.adjust,size = Count))+
  scale_color_gradient(low ="#FF0000",high = "#00A08A")+
  xlab("Gene Ratio")+
  theme_bw()+
  guides(color = guide_colorbar(reverse = TRUE))
dev.off()
write.csv(test,file = "feature.GO.BP.BH.csv")
write.csv(test,file = "KEGG.csv")

results<-read.table("feature-GOBP-filter.txt",sep = "\t",header = T,row.names = 1)
results<-results[order(results$fold,decreasing = F),]
results$Description<-factor(results$Description,levels = results$Description)

jpeg(file="feature-GOBP-filter.jpeg",width =10,height =6,units = "in", res = 1000)
ggplot(results,aes(x=fold,y=Description,fill=pvalue))+geom_bar(stat="identity")+
  scale_fill_continuous(low="#B40F20" ,high="#5BBCD6")+
  theme_bw()+theme(        axis.text.x = element_text(size = 15),
                           axis.ticks.x = element_blank(),
                           axis.text.y = element_text(size = 15))
dev.off()
#################feature genes heatmap#############
setwd("D:/Bioinfrolf/Bioinfrolf/GSEA/")
Genes<-read.table("D:/Bioinfrolf/Bioinfrolf/model/impMatrix.txt",header = T)
genes=Genes$Feature#fNamesAd#

load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
Mono_Macr<-subset(x = Mono_Macr,subset = ann!="Macr-CCL5")
Mono_Macr@meta.data$ann<-Mono_Macr@meta.data$ann[,drop=T]
clusters<-Mono_Macr@meta.data$ann%>%as.factor()%>%levels()
data<-Mono_Macr@assays$RNA@data
data<-data[genes,]%>%as.matrix()

ann_row=data.frame(CellType=Mono_Macr@meta.data$ann,
                   Tissue=Mono_Macr@meta.data$tissue_sub,
                   row.names = rownames(Mono_Macr@meta.data))
head(ann_row)
table(ann_row)
d1<-data[,rownames(ann_row)[ann_row$CellType==clusters[1]]]
d2<-data[,rownames(ann_row)[ann_row$CellType==clusters[2]]]
d3<-data[,rownames(ann_row)[ann_row$CellType==clusters[3]]]
d4<-data[,rownames(ann_row)[ann_row$CellType==clusters[4]]]
d5<-data[,rownames(ann_row)[ann_row$CellType==clusters[5]]]
d6<-data[,rownames(ann_row)[ann_row$CellType==clusters[6]]]
d7<-data[,rownames(ann_row)[ann_row$CellType==clusters[7]]]
# d8<-data[,rownames(ann_row)[ann_row$CellType==clusters[8]]]
core<-data[,rownames(ann_row)[ann_row$Tissue=="TumorCore"]]
edge<-data[,rownames(ann_row)[ann_row$Tissue=="TumorEdge"]]
normal<-data[,rownames(ann_row)[ann_row$Tissue=="Normal"]]

data<-cbind(d1,d2,d3,d4,d5,d6,d7)
tdata<-cbind(normal,edge,core)

d9<-apply(normal,1, mean)
d10<-apply(edge,1, mean)
d11<-apply(core,1, mean)
d1<-apply(d1, 1, mean)
d2<-apply(d2, 1, mean)
d3<-apply(d3, 1, mean)
d4<-apply(d4, 1, mean)
d5<-apply(d5, 1, mean)
d6<-apply(d6, 1, mean)
d7<-apply(d7, 1, mean)
# d8<-apply(d8, 1, mean)

ors = list(
  CellType = c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
               "Macr-BPL5"="#B40F20","Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
               "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A"),
  Tissue=c("Adjacent Normal"="#00A08A", "Peripheral Tumor"="#F98400","Core Tumor"="#FF0000")
)

df<-cbind(d1,d2,d3,d4,d5,d6,d7)
colnames(df)<-clusters
df<-as.data.frame(df)
df$gene_name<-rownames(df)
rownames(df)<-1:nrow(df)
cirdf <- melt(df)

#################circle tissue heatmap#######
library(gplots)
library(ggdendro)
library(ggh4x)

fdf<-cbind(d9,d10,d11)
colnames(fdf)<-c("Adjacent Normal","Peripheral Tumor", "Core Tumor")
head(fdf)
class(df)


colorRampPalette(c('#FF6767',"white"))(10)
colorRampPalette(c("#00A08A","#F98400","#FF0000"))(11)
color=colorRampPalette(c("#313695","white","#A50026"))(100)

h<-pheatmap::pheatmap(fdf,scale = "row",
                      show_colnames =T,
                      show_rownames = T,
                      cluster_row = T, cluster_col = T,
                      border_color =NA,
                      #kmeans_k = 3,
                      #clustering_method = "ward.D2",
                      color = colorRampPalette(
                        rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                      angle_col=45,fontsize=10,fontsize_col=10,
                      main ="Feature Expression")

jpeg(file="tissue.jpeg",width=9.5,height=8,units = "in", res = 2000)
h
dev.off()

name<-rownames(fdf)[h$tree_row$order]

####
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
mat<-scale_rows(fdf)
rownames(mat)<-rownames(fdf)
colnames(mat)<-colnames(fdf)
mat<-mat[name,]

mat<-as.data.frame(mat)
mat$gene_name<-rownames(mat)
rownames(mat)<-1:nrow(mat)
cirmat<-melt(mat)
head(cirmat)

resn <- cirmat %>% filter(variable == 'Adjacent Normal')
rese <- cirmat %>% filter(variable == 'Peripheral Tumor')
resc <- cirmat %>% filter(variable == 'Core Tumor')

resn$ang <- seq(from = (360/nrow(resn)) / 1.5,
                to = (1.5* (360/nrow(resn))) - 360,
                length.out = nrow(resn)) + 80

resn$hjust <- 0
resn$hjust[which(resn$ang < -90)] <- 1
resn$ang[which(resn$ang < -90)] <- (180+resn$ang)[which(resn$ang < -90)]
resn$labels<-NA
resn$labels[as.numeric(rownames(resn))%%5==1]<-resn$gene_name[as.numeric(rownames(resn))%%10==1]

cirmat$var <- rep(1:nrow(resn),3)
head(cirmat)
dim(cirmat)

g_anno2 <- data.frame(importance = Genes$Gain,Genes=Genes$Feature)
rownames(g_anno2)<-Genes$Feature
g_anno2<-g_anno2[resn$gene_name,]

g_anno1 <- data.frame(Genes=Genes$Feature,
                      `Highly Expressed`=ifelse(Genes$Feature%in%G0,"Adjacent Normal",
                                                ifelse(Genes$Feature%in%G1,"Peripheral Tumor","Core Tumor")))


# g_anno1 <- data.frame(type = c(rep('Normal',200),rep('Core Tumor High Expression',100),
#                                rep('Peripheral Tumor Highly Expressed',100),rep('Normal High Expression',45)))

s_anno <- data.frame(value = c("Adjacent Normal", "Peripheral Tumor","Core Tumor"))
s_anno$value <-factor(s_anno$value,levels = c("Adjacent Normal", "Peripheral Tumor","Core Tumor"))
cirmat$variable<-factor(cirmat$variable,levels = c("Adjacent Normal", "Peripheral Tumor","Core Tumor"))

scale_fill_gradient(low = "#FFEEEE",high = '#FF6767')
b<-seq(1,445,5)
p1 <- ggplot() +
  geom_bar(data = g_anno2,stat = 'identity',
           aes(x = 1:nrow(resn),y = -0.5,fill = importance),
           width = 1,
           color = NA) +
  scale_fill_gradient(low = "white",high = '#B40F20')+
  scale_x_discrete(breaks = b)+
  new_scale("fill") +
  geom_bar(data = g_anno1,stat = 'identity',
           aes(x = 1:nrow(resn),y = -0.25,fill = Highly.Expressed),
           width = 1,
           color = NA) +
  scale_fill_manual(name = 'Highly Expressed',
                    values =
                      c('Adjacent Normal'="#00A08A",
                        'Peripheral Tumor'="#F97500",
                        'Core Tumor'="#FC3400"))+
  new_scale("fill") +
  geom_bar(data = s_anno,stat = 'identity',
           aes(x = -2,y = 1,fill = value),
           width = 4,
           color = NA) +
  scale_fill_manual(name = 'Tissue',
                    values = c("Adjacent Normal"="#00A08A", 
                               "Peripheral Tumor"="#F98400",
                               "Core Tumor"="#FF0000"))

p3 <- p1 + 
  new_scale("fill") +
  geom_tile(data = cirmat[which(cirmat$variable == 'Adjacent Normal'),],
            aes(x = 1:nrow(resn),y = 2.5,fill = value),
            color = 'white') +
  geom_tile(data = cirmat[which(cirmat$variable == 'Peripheral Tumor'),],
            aes(x = 1:nrow(resn),y = 1.5,fill = value),
            color = 'white') +
  geom_tile(data = cirmat[which(cirmat$variable == 'Core Tumor'),],
            aes(x = 1:nrow(resn),y = 0.5,fill = value),
            color = 'white') +
  scale_fill_gradient2(midpoint = 0,
                       low = "#5BBCD6",
                       mid = "white",
                       high = "#B40F20") +
  ylim(-3,5)+xlim(-20,nrow(resn)+5)

jpeg(file="tissue_3.jpeg",width=15,height=10,units = "in", res = 1000)
p3 + coord_polar(theta = 'x') +theme_void() +
  geom_text(data = resn,
            aes(x = as.numeric(rownames(resn)),
                y = 3.2,
                label = labels, angle = ang, hjust = hjust,family="serif"),#
            size = 4)#+NoLegend()
dev.off()

###
head(mat)
dac <- mat[,-4]
rownames(dac) <- mat$gene_name

yclust <- hclust(dist(dac))
xclust <- hclust(dist(t(dac)))

jpeg(file="tissue_1.jpeg",width=4,height=15,units = "in", res = 2000)
ggplot(cirmat,aes(x = variable,y = gene_name)) +
  geom_tile(aes(fill = value),color = 'white') +
  scale_fill_gradient2(midpoint = 0,
                       low = "#5BBCD6",
                       mid = "white",
                       high = "#B40F20") +
  scale_y_dendrogram(hclust = yclust) +
  scale_x_dendrogram(hclust = xclust,position = 'top') +
  theme_classic() +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 2)) +
  xlab('') + ylab('')
dev.off()

# set.seed(1231)
# ind <- sample(2, ncol(data), replace = TRUE, prob = c(0.1, 0.9))
# test <- data[,ind==1] #the test data set


