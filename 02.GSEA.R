library("org.Hs.eg.db")
library(enrichplot)
library(msigdbr)
library("RColorBrewer")
library(wesanderson)
library(clusterProfiler)
library(scales)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(ggplot2)
library(stringr)
library(reshape2)
library(rlang)
library(Seurat)
library(fgsea)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(ggrepel)
library(ggnewscale)
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/03.MNPs.Rdata")
passway.1<-read.table("D:/Bioinfrolf/Bioinfrolf/data/GSVA_3.txt")
passway.p<-read.table("D:/Bioinfrolf/Bioinfrolf/data/GSVA_p.txt",header = F)
features = NULL
ribo <- features %||% grep(pattern = "^RP[SL]", x = rownames(sce@assays$RNA@counts),value = TRUE)
mt<- features %||% grep(pattern = "^MT-", x = rownames(sce@assays$RNA@counts),value = TRUE)
counts <- GetAssayData(object = sce, slot = "counts")
filtered_counts <- counts[!rownames(counts)%in%c(ribo,mt),]
pbmc <- CreateSeuratObject(filtered_counts, meta.data = sce@meta.data)%>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData(verbose = FALSE)

clusters<-sce@meta.data$ann%>%as.factor()%>%levels()
msigdbr_species() 
hs_msigdbr <- msigdbr(species="Homo sapiens")
hsGO <- msigdbr(species="Homo sapiens",category="C5")
hsGOBP <- msigdbr(species="Homo sapiens",category="C5",subcategory = "BP")
head(hsGO)%>%as.data.frame()
unique(hsGOBP$gs_exact_source)%>%head()
unique(hsGOBP$gs_subcat)

library(presto)
pbmc.genes <- wilcoxauc(pbmc, assay = "data",'ann')
head(pbmc.genes)
dplyr::count(pbmc.genes, group)
#################################
fgsea_sets<- hsGOBP %>% split(x = .$gene_symbol, f = .$gs_name)
runGSEA<-function(pbmc.genes,cluster){
  cluster.genes<- pbmc.genes %>% dplyr::filter(group == cluster) %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
  ranks<- deframe(cluster.genes)
  head(ranks)
  library(msigdbr)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
  fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
  fgseaResTidy$clusters<-rep(cluster,nrow(fgseaResTidy))
  return(fgseaResTidy)
}

res1<-runGSEA(pbmc.genes,clusters[1])#cluster<-"Macr-APOE"
res2<-runGSEA(pbmc.genes,clusters[2])#"Macr-RPs"
res3<-runGSEA(pbmc.genes,clusters[3])#"Mono-VCAN"
res4<-runGSEA(pbmc.genes,clusters[4])#"Mono-CD16"
res5<-runGSEA(pbmc.genes,clusters[5])#"Macr-CXCL10"
res6<-runGSEA(pbmc.genes,clusters[6])
res7<-runGSEA(pbmc.genes,clusters[7])#"Macr-VCAM1"

save(hsGO,hsGOBP,pbmc.genes,clusters,file = "03.MNPs_GSVA.Rdata")

filterRES<-function(res){
  res<-res[res$padj<0.05,]
  res<-res[order(abs(res$NES),decreasing = T),]
  res<-as.data.frame(res)
  rownames(res)<-res$pathway
  res$pathway<-lapply(rownames(res),function(x){
    paste(unlist(strsplit(x,split = "_"))[-1],
          collapse = " ")%>% str_to_title()
  })%>%unlist()
  a<-gsub("In","in",res$pathway)
  a<-gsub("Of","of",a)
  res$pathway<-gsub("To","to",a)
  return(res)
}
res1f<-filterRES(res1)
res2f<-filterRES(res2)
res3f<-filterRES(res3)
res4f<-filterRES(res4)
res5f<-filterRES(res5)
res6f<-filterRES(res6)
res7f<-filterRES(res7)

g<-intersect(res1f$pathway,res2f$pathway)
g<-intersect(g,res3f$pathway)
g<-intersect(g,res4f$pathway)
g<-intersect(g,res5f$pathway)
g<-intersect(g,res6f$pathway)
g<-intersect(g,res7f$pathway)
g<-passway.p$V2
df<-rbind(res1f[g,],res2f[g,],res3f[g,],res4f[g,],res5f[g,],res6f[g,],res7f[g,])
write.csv(res1f[g,1:4],file = "GSVA_p.csv",quote = F)
write.csv(res1f[,-8],file = "GSVA_Macr-APOE.csv")
write.csv(res2f[,-8],file = "GSVA_Macr-RPs.csv")
write.csv(res3f[,-8],file = "GSVA_Mono-VCAN.csv")
write.csv(res4f[,-8],file = "GSVA_Mono-CD16.csv")
write.csv(res5f[,-8],file = "GSVA_Macr-STAT1.csv")
write.csv(res6f[,-8],file = "GSVA_Macr-SPP1.csv")
write.csv(res7f[,-8],file = "GSVA_Macr-VCAM1.csv")

{
  intersect(g,passway$V1)
  
  filterRES<-function(res){
    res<-res[order(abs(res$NES),decreasing = T),]
    res<-as.data.frame(res)
    rownames(res)<-res$pathway
    res$pathway<-lapply(rownames(res),function(x){
      paste(unlist(strsplit(x,split = "_"))[-1],
            collapse = " ")%>% str_to_title()
    })%>%unlist()
    a<-gsub("In","in",res$pathway)
    a<-gsub("Of","of",a)
    res$pathway<-gsub("To","to",a)
    return(res)
  }
  g<-passway.1$V1%>%unique()
  df<-rbind(res1f[g,],res2f[g,],res3f[g,],res4f[g,],res5f[g,],res6f[g,],res7f[g,])
  data<-df[,c(1,5,9)]%>%as.data.frame()
  library(reshape2)
  data<-melt(data,id.vars = c("clusters", "pathway"))
  data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
  rownames(data)<-data$pathway
  data<-data[,-1]
  data<-as.matrix(data)
  matBP<-data
  
  data<-df[,c(1,3,9)]%>%as.data.frame()
  library(reshape2)
  data<-melt(data,id.vars = c("clusters", "pathway"))
  data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
  rownames(data)<-data$pathway
  data<-data[,-1]
  data<-as.matrix(data)
  matpval<-data
}

data<-df[,c(1,5,9)]%>%as.data.frame()
library(reshape2)
data<-melt(data,id.vars = c("clusters", "pathway"))
data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
rownames(data)<-data$pathway
data<-data[,-1]
data<-as.matrix(data)
head(data)
rownames(data)

matBP<-data
rownames(matBP)<-lapply(rownames(matBP),function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()
head(matBP)

a<-gsub("In","in",rownames(matBP))
a<-gsub("Of","of",a)
rownames(matBP)<-gsub("To","to",a)

data<-df[,c(1,2,9)]%>%as.data.frame()
library(reshape2)
data<-melt(data,id.vars = c("clusters", "pathway"))
data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
rownames(data)<-data$pathway
data<-data[,-1]
data<-as.matrix(data)
head(data)
rownames(data)

matpval<-data
rownames(matpval)<-lapply(rownames(matpval),function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()
head(matpval)

a<-gsub("In","in",rownames(matpval))
a<-gsub("Of","of",a)
rownames(matpval)<-gsub("To","to",a)
library(pheatmap)

max(df$pval)
min(df$NES)
matpval[1:5,1:5]
matplot<-matrix(ifelse(matpval >= 0.05, 0, matBP),nrow(matBP))
matplot1<-matrix(ifelse(matpval >= 0.05, NA, matBP),nrow(matBP))
rownames(matplot)<-rownames(matBP)
colnames(matplot)<-colnames(matBP)
rownames(matplot1)<-rownames(matBP)
colnames(matplot1)<-colnames(matBP)
matplot[1:5,1:5]
p<-matpval[which(apply(matplot, 1,sd)!=0),]
p1<-matBP[which(apply(matplot, 1,sd)!=0),]
p2<-p1[order(apply(p1, 1,sd),decreasing = T),]
head(p2)
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(1,0.05,0), c("light grey","white","#FF0000"))
col_fun1= colorRamp2(c(min(df$NES),0,max(df$NES)), c("#5BBCD6","white","#B40F20"))
int<-abs(max(df$NES)-min(df$NES))/100
bk <- c(seq(min(df$NES),-int,by=int),seq(0,max(df$NES),by=int))

cols<-c(colorRampPalette(colors = c("#5BBCD6","white"))(which(bk==0)),
        colorRampPalette(colors = c("white","#B40F20"))(100-which(bk==0)))
#colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100)
cn<-c("Macr-VCAM1","Macr-APOE", "Macr-SPP1" , "Macr-STAT1", "Macr-RPs", "Mono-CD16" , "Mono-VCAN" )
matBP<-matBP[unique(df$pathway),cn]
matpval<-matpval[unique(df$pathway),cn]
colnames(matBP)

ann_row=data.frame(PasswayType=rep(c("Inflammation","Antigen Presentation","Angiogenesis",
                                  "Remodeling Cell Adhesion and Matrix","Migration","Oxidative Stress"),
                                c(11,6,3,5,8,5)),
                   row.names = rownames(matBP))
jpeg(file="heatmap_GSVA_all-4.jpeg",width=12,height=9,units = "in", res = 1000)
Heatmap(matBP, 
        row_split = ann_row$PasswayType, 
        #column_split = ann_row$CellType,
        cluster_rows = FALSE,
        cluster_columns = T,
        show_column_names = T,
        show_row_names =T,
        col = col_fun1
)
dev.off()
jpeg(file="heatmap_GSVA_all-4.jpeg",width=10,height=9,units = "in", res = 1000)
pheatmap::pheatmap(matBP,
                   number_color="black",
                   display_numbers = matrix(ifelse(matpval< 0.001, "***", ifelse(
                     matpval < 0.01,"**",ifelse(
                       matpval< 0.05,"*",""
                     )
                   )),nrow(matpval)),
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = F, cluster_col = T,
                   border_color =NA,
                   angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,fontsize_number=15,
                   color = cols,
                   main ="HCC MNPs GSVA")
dev.off()
write.csv(df[1:96,c(1,3,5,9)],file = "GSVA_1.csv",quote = F)
####################MNPs####################
setwd("D:/Bioinfrolf/Bioinfrolf/GSEA/")
macrolist<-read.delim("macrophage.txt",header = F)
monolist<-read.delim("monocyte.txt",header = F)
list<-rbind(macrolist,monolist)
list<-list[!duplicated(list$V1),]
rownames(list)<-list$V1
m_df<- hsGO[hsGO$gs_exact_source%in%list$V1,]
head(m_df)
library(fgsea)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(fgsea_sets)

runGSEA<-function(pbmc.genes,cluster){
  cluster.genes<- pbmc.genes %>% dplyr::filter(group == cluster) %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
  ranks<- deframe(cluster.genes)
  head(ranks)
  library(msigdbr)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
  fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
  fgseaResTidy$clusters<-rep(cluster,nrow(fgseaResTidy))
  return(fgseaResTidy)
}

res1<-runGSEA(pbmc.genes,clusters[1])#cluster<-"Macr-APOE"
res2<-runGSEA(pbmc.genes,clusters[2])#"Macr-RPs"
res3<-runGSEA(pbmc.genes,clusters[3])#"Mono-VCAN"
res4<-runGSEA(pbmc.genes,clusters[4])#"Mono-CD16"
res5<-runGSEA(pbmc.genes,clusters[5])#"Macr-CXCL10"
res6<-runGSEA(pbmc.genes,clusters[6])
res7<-runGSEA(pbmc.genes,clusters[7])#"Macr-VCAM1"

df<-rbind(res1,res2,res3,res4,res5,res6,res7)
# write.csv(df,file = "MNPs_GSVAres.csv",quote = F)
# save(df,file = "MNPs_GSVAres.Rdata")

data<-df[,c(1,5,9)]%>%as.data.frame()
library(reshape2)
data<-melt(data,id.vars = c("clusters", "pathway"))
data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
rownames(data)<-data$pathway
data<-data[,-1]
data<-as.matrix(data)
head(data)
rownames(data)

matBP<-data[1:59,]
rownames(matBP)<-lapply(rownames(matBP),function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()
head(matBP)

matBP<-matBP[-grep(pattern ="Osteoclast",rownames(matBP)),]
matBP<-matBP[-grep(pattern ="Marginal Zone B",rownames(matBP)),]
a<-gsub("In","in",rownames(matBP))
a<-gsub("Of","of",a)
rownames(matBP)<-gsub("To","to",a)

data<-df[,c(1,2,9)]%>%as.data.frame()
library(reshape2)
data<-melt(data,id.vars = c("clusters", "pathway"))
data<-dcast(data, pathway ~ clusters, na.rm = TRUE)
rownames(data)<-data$pathway
data<-data[,-1]
data<-as.matrix(data)
head(data)
rownames(data)

matpval<-data[1:59,]
rownames(matpval)<-lapply(rownames(matpval),function(x){
  paste(unlist(strsplit(x,split = "_"))[-1],
        collapse = " ")%>% str_to_title()
})%>%unlist()
head(matpval)

matpval<-matpval[-grep(pattern ="Osteoclast",rownames(matpval)),]
matpval<-matpval[-grep(pattern ="Marginal Zone B",rownames(matpval)),]
a<-gsub("In","in",rownames(matpval))
a<-gsub("Of","of",a)
rownames(matpval)<-gsub("To","to",a)
library(pheatmap)

max(df$pval)
min(df$NES)
matpval[1:5,1:5]
matplot<-matrix(ifelse(matpval >= 0.05, 0, matBP),nrow(matBP))
matplot1<-matrix(ifelse(matpval >= 0.05, NA, matBP),nrow(matBP))
rownames(matplot)<-rownames(matBP)
colnames(matplot)<-colnames(matBP)
rownames(matplot1)<-rownames(matBP)
colnames(matplot1)<-colnames(matBP)
matplot[1:5,1:5]
p<-matpval[which(apply(matplot, 1,sd)!=0),]
p1<-matBP[which(apply(matplot, 1,sd)!=0),]
p2<-p1[order(apply(p1, 1,sd),decreasing = T),]
head(p2)
library(circlize)
library(ComplexHeatmap)
col_fun = colorRamp2(c(1,0.05,0), c("light grey","#FF0000"))#"white",

jpeg(file="heatmap_GSVA_MNPs.jpeg",width=10,height=8,units = "in", res = 1000)
pheatmap::pheatmap(p1,
                   #scale = "row",
                   number_color="white",
                   display_numbers = matrix(ifelse(p< 0.001, "***", ifelse(
                     p < 0.01,"**",ifelse(
                       p< 0.05,"*",""
                     )
                   )),nrow(p1)),
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = T, cluster_col = T,
                   border_color =NA,
                   #fontface="bold",
                   angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,fontsize_number=15,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                             "RdYlBu")))(100),
                   main ="HCC MNPs GSVA")
dev.off()

########################
results<-read.table("feature-GOBP-filter.txt",sep = "\t",header = T,row.names = 1)
results<-results[order(results$p.adjust,decreasing = T),]
results1<-results[-grep(pattern ="migration",results$Description),]
results2<-results[grep(pattern ="migration",results$Description),]
results<-results[results$p.adjust<=0.05,]
results$Description<-factor(results$Description,levels = results$Description)
#results$pvalue[5]=0.05
jpeg(file="feature-GOBP-filter-2.jpeg",width =6,height =4,units = "in", res = 1000)
ggplot(results2,aes(x=fold,y=Description,fill=p.adjust))+geom_bar(stat="identity")+
  scale_fill_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
  theme_bw()+theme(        axis.text.x = element_text(size = 12),
                           axis.text.y = element_text(size = 12))+
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 15, color = 'black'),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_text(size=12), 
        legend.key.size=unit(0.2, "inches")
  )+labs(title ="",x = 'NES',y = "")
dev.off()
