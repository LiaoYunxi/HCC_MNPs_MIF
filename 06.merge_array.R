rm(list = ls())
#BiocManager::install("BatchQC")
library(GEOquery)
library(limma)
library(umap)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(BatchQC)
###############illumina#################
{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\GSE121248-HCC-array")
  exprMatrix.1 = read.table(file = "GSE121248_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.1<-data.table::fread("meta.txt",header = F)
  meta.1<-meta.1[order(meta.1$V2,decreasing = T),]
  
  setwd("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/GSE41804-HCC-array")
  exprMatrix.2 = read.table(file = "GSE41804_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.2<-read.table("meta.txt")
  meta.2<-meta.2[order(meta.2$V3,decreasing = F),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB8_mRNA")
  exprMatrix.3 = read.table(file = "HCCDB8_GSE9843_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.3<-data.table::fread(file = "HCCDB8.sample.txt",header = F)
  meta.3<-t(meta.3)%>%as.data.frame()
  colnames(meta.3)<-meta.3[1,]
  meta.3<-meta.3[-1,]
  
  index=intersect(intersect(rownames(exprMatrix.1),rownames(exprMatrix.2)),rownames(exprMatrix.3))
  exprMatrix.1 =exprMatrix.1[index,meta.1$V1]
  exprMatrix.2=exprMatrix.2[index,meta.2$V1]
  exprMatrix.3=exprMatrix.3[index,meta.3$GEO_ID]
  
  
  tumor.1<-exprMatrix.1[,1:70]
  normal.1<-exprMatrix.1[,71:ncol(exprMatrix.1)]
  which(meta.2$V3=="non-tumor1")
  tumor.2<-exprMatrix.2[,1:20]
  normal.2<-exprMatrix.2[,21:ncol(exprMatrix.2)]
  tumor.3<-exprMatrix.3
  
  exprMatrix<-cbind(exprMatrix.1,exprMatrix.2,exprMatrix.3)
  tumor<-cbind(tumor.1,tumor.2,tumor.3)
  normal<-cbind(normal.1,normal.2)
}
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
save(tumor,normal,exprMatrix,data,ex_b_limma,ex_b_sva,file = "illumina1.Rdata")
load("illumina1.Rdata")

{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB4_mRNA")
  exprMatrix.1 = read.table(file = "HCCDB4_GSE36376_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  n<-which(is.na(exprMatrix.1[1,]))
  exprMatrix.1<-exprMatrix.1[,-n]
  # attr(na.omit(exprMatrix.1),"na.action")
  # exprMatrix.1<-na.omit(exprMatrix.1)
  meta.1<-data.table::fread(file = "HCCDB4.sample.txt",header = F)
  meta.1<-t(meta.1)%>%as.data.frame()
  colnames(meta.1)<-meta.1[1,]
  meta.1<-meta.1[-1,]
  rownames(meta.1)<-meta.1$GEO_ID
  meta.1<-meta.1[colnames(exprMatrix.1),]
  meta.1<-meta.1[order(meta.1$TYPE,decreasing = T),]
  
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB14_mRNA")
  exprMatrix.2 = read.table(file = "HCCDB14_GSE43619_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.2<-data.table::fread(file = "HCCDB14.sample.txt",header = F)
  meta.2<-t(meta.2)%>%as.data.frame()
  colnames(meta.2)<-meta.2[1,]
  meta.2<-meta.2[-1,]
  rownames(meta.2)<-meta.2$GEO_ID
  meta.2<-meta.2[order(meta.2$TYPE,decreasing = T),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB17_mRNA")
  exprMatrix.3 = read.table(file = "HCCDB17_GSE76427_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.3<-data.table::fread(file = "HCCDB17.sample.txt",header = F)
  meta.3<-t(meta.3)%>%as.data.frame()
  colnames(meta.3)<-meta.3[1,]
  meta.3<-meta.3[-1,]
  rownames(meta.3)<-meta.3$GEO_ID
  meta.3<-meta.3[order(meta.3$TYPE,decreasing = T),]
  max(which(meta.3$TYPE=="HCC"))
  

  index=intersect(intersect(rownames(exprMatrix.1),rownames(exprMatrix.2)),rownames(exprMatrix.3))
  exprMatrix.1 =exprMatrix.1[index,meta.1$GEO_ID]
  exprMatrix.2=exprMatrix.2[index,meta.2$GEO_ID]
  exprMatrix.3=exprMatrix.3[index,meta.3$GEO_ID]
  
  
  tumor.1<-exprMatrix.1[,1:max(which(meta.1$TYPE=="HCC"))]
  normal.1<-exprMatrix.1[,(max(which(meta.1$TYPE=="HCC"))+1):ncol(exprMatrix.1)]
  tumor.2<-exprMatrix.2
  tumor.3<-exprMatrix.3[,1:max(which(meta.3$TYPE=="HCC"))]
  normal.3<-exprMatrix.3[,(max(which(meta.3$TYPE=="HCC"))+1):ncol(exprMatrix.3)]
  
  exprMatrix<-cbind(exprMatrix.1,exprMatrix.2,exprMatrix.3)
  tumor<-cbind(tumor.1,tumor.2,tumor.3)
  normal<-cbind(normal.1,normal.3)
  data<-cbind(normal,tumor)
}
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
save(tumor,normal,exprMatrix,data,ex_b_sva,file = "illumina2.Rdata")
load("illumina2.Rdata")

{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB9_mRNA")
  exprMatrix.1 = read.table(file = "HCCDB9_GSE19977_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.1<-data.table::fread(file = "HCCDB9.sample.txt",header = F)
  meta.1<-t(meta.1)%>%as.data.frame()
  colnames(meta.1)<-meta.1[1,]
  meta.1<-meta.1[-1,]
  meta.1<-meta.1[order(meta.1$TYPE,decreasing = T),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB11_mRNA")
  exprMatrix.2 = read.table(file = "HCCDB11_GSE46444_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.2<-data.table::fread(file = "HCCDB11.sample.txt",header = F)
  meta.2<-t(meta.2)%>%as.data.frame()
  colnames(meta.2)<-meta.2[1,]
  meta.2<-meta.2[-1,]
  rownames(meta.2)<-meta.2$GEO_ID
  meta.2<-meta.2[order(meta.2$TYPE,decreasing = T),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
  load("illumina1.Rdata")
  load("illumina2.Rdata")
  index=intersect(intersect(rownames(exprMatrix.1),rownames(ex_b_sva)),rownames(ex_b_limma))
  exprMatrix.1 =exprMatrix.1[index,meta.1$GEO_ID]
  exprMatrix.2=exprMatrix.2[index,meta.2$GEO_ID]
  exprMatrix.3=ex_b_limma[index,]%>%as.data.frame()
  exprMatrix.4=ex_b_sva[index,]%>%as.data.frame()
  
  tumor.1<-exprMatrix.1
  tumor.2<-exprMatrix.2[,1:max(which(meta.2$TYPE=="HCC"))]
  normal.2<-exprMatrix.2[,(max(which(meta.2$TYPE=="HCC"))+1):ncol(exprMatrix.2)]
  tumor.3<-exprMatrix.3[,58:238]
  dim(ex_b_limma)
  normal.3<-exprMatrix.3[,1:57]
  tumor.4<-exprMatrix.4[,246:684]
  dim(ex_b_sva)
  normal.4<-exprMatrix.4[,1:245]
  
  exprMatrix<-cbind(exprMatrix.1,exprMatrix.2,exprMatrix.3,exprMatrix.4)
  tumor<-cbind(tumor.1,tumor.2,tumor.3,tumor.4)
  normal<-cbind(normal.2,normal.3,normal.4)
  data<-cbind(normal,tumor)
}

{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
  load("illumina1.Rdata")
  load("illumina2.Rdata")
  index=intersect(rownames(ex_b_limma),rownames(ex_b_sva))
  exprMatrix.3=ex_b_limma[index,]%>%as.data.frame()
  exprMatrix.4=ex_b_sva[index,]%>%as.data.frame()
  
  tumor.3<-exprMatrix.3[,58:238]
  dim(ex_b_limma)
  normal.3<-exprMatrix.3[,1:57]
  tumor.4<-exprMatrix.4[,246:684]
  dim(ex_b_sva)
  normal.4<-exprMatrix.4[,1:245]
  
  exprMatrix<-cbind(exprMatrix.3,exprMatrix.4)
  tumor<-cbind(tumor.3,tumor.4)
  normal<-cbind(normal.3,normal.4)
  data<-cbind(normal,tumor)
}
rt=cbind(ID=rownames(ex_b_sva),as.data.frame(ex_b_sva))
write.table(rt,file = "illumina_symbol_tpm.txt",quote = F,row.names = F)
###############remove Batch Effect###################
ex<-exprMatrix
{
  batch1=c(rep(c("GSE121248-GPL570","GSE41804-GPL570"),c(37,20)),
           rep(c("GSE121248-GPL570","GSE41804-GPL570","GSE9843-GPL570"),c(70,20,91)))
  tissue1<-rep(c("Normal","Tumor"),c(57,181))
  data<-cbind(normal,tumor)
}

{
  batch=c(rep(c("GSE36376","GSE76427"),c(193,52)),
          rep(c("GSE36376","GSE43619","GSE76427"),c(236,88,115)))
  
  tissue<-rep(c("Normal","Tumor"),c(245,439))
}

{
  batch=c(rep(c("GSE46444","GPL570","GPL10558"),c(48,57,245)),
          rep(c("GSE19977","GSE46444","GPL570","GPL10558"),c(164,88,181,439)))
  
  tissue<-rep(c("Normal","Tumor"),c(350,872))
}

{
  batch=c(rep(c("GPL570","GPL10558"),c(57,245)),
          rep(c("GPL570","GPL10558"),c(181,439)))
  
  tissue<-rep(c("Normal","Tumor"),c(57+245,181+439))
}
# which(is.na(data))
# data<-na.omit(data)
design <- model.matrix(~0+as.factor(tissue))
colnames(design)<-c("Normal","Tumor")
rownames(design)<-colnames(data)

ex_b_limma <- removeBatchEffect(data,
                                batch = batch,design = design)
library(sva)
design <- model.matrix(~+as.factor(tissue1))
ex_b_sva = ComBat(dat=as.matrix(data), 
                  batch=batch1,mod = design)
illumina1<-ex_b_sva
{
  library(DESeq2)
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                              colData = colData,  design= ~ batch)
  
  dds <- DESeq(ddsFullCountTable)
  
  colData <- data.frame(row.names=colnames(data), batch)
  dds <- estimateSizeFactorsForMatrix(data) 
  # idx <- rowMeans(data) > 1
  # data <- data[idx, ]
  mod <- model.matrix(~ factor(batch), colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 2)
}

e<-expr
group_list<-batch
colnames(e)<-group_list
library(factoextra)
e<-t(e)
e<-as.data.frame(e)
res.pca <- prcomp(e, scale = TRUE)

#jpeg(file="feature_pca.jpeg",width =6,height =5,units = "in", res = 2000)
fviz_pca_ind(res.pca,
             col.ind = group_list, # ??ɫ??Ӧgroup??Ϣ
             #palette =wes_palette(n=2, name="BottleRocket2"),
             #addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Group",## Legend????
             repel = TRUE
)
#dev.off()

ex<-exprMatrix
{
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) }
  
  par(mar=c(7,4,2,1))
  title <- paste ("illumina", "/", "GPL570", sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
  dev.off()
  
  # expression value distribution plot
  par(mar=c(4,4,2,1))
  plotDensities(ex, main=title, legend=F)
  
  # mean-variance trend
  ex <- na.omit(ex) # eliminate rows with NAs
  plotSA(lmFit(ex), main="Mean variance trend, GPL570")
  
  # UMAP plot (multi-dimensional scaling)
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
  df<-as.data.frame(ump$layout)
  ggplot(df,aes(x=V1,y=V2))+geom_point(aes(color = batch),size=2)+
    theme_bw()+theme(panel.grid = element_blank(),
                     axis.line = element_line(size = 0.5),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_text(face = 'bold',size = 5),
    )+xlab('UMAP1') + ylab('UMAP2')

}

###############affy#############
{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB1_mRNA")
  exprMatrix.1 = read.table(file = "HCCDB1_GSE22058_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.1<-data.table::fread(file = "HCCDB1.sample.txt",header = F)
  meta.1<-t(meta.1)%>%as.data.frame()
  colnames(meta.1)<-meta.1[1,]
  meta.1<-meta.1[-1,]
  rownames(meta.1)<-meta.1$GPL6793_GEO_ID
  meta.1<-meta.1[order(meta.1$TYPE,decreasing = T),]
  
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB3_mRNA")
  exprMatrix.2 = read.table(file = "HCCDB3_GSE25097_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.2<-data.table::fread(file = "HCCDB3.sample.txt",header = F)
  meta.2<-t(meta.2)%>%as.data.frame()
  colnames(meta.2)<-meta.2[1,]
  meta.2<-meta.2[-1,]
  rownames(meta.2)<-meta.2$GEO_ID
  meta.2<-meta.2[order(meta.2$TYPE,decreasing = T),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB6_mRNA")
  exprMatrix.3 = read.table(file = "HCCDB6_GSE14520-GPL3921_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.3<-data.table::fread(file = "HCCDB6.sample.txt",header = F)
  meta.3<-t(meta.3)%>%as.data.frame()
  colnames(meta.3)<-meta.3[1,]
  meta.3<-meta.3[-1,]
  rownames(meta.3)<-meta.3$GEO_ID
  meta.3<-meta.3[order(meta.3$TYPE,decreasing = T),]
  max(which(meta.3$TYPE=="HCC"))
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB13_mRNA")
  exprMatrix.4 = read.table(file = "HCCDB13_GSE63898_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.4<-data.table::fread(file = "HCCDB13.sample.txt",header = F)
  meta.4<-t(meta.4)%>%as.data.frame()
  colnames(meta.4)<-meta.4[1,]
  meta.4<-meta.4[-1,]
  rownames(meta.4)<-meta.4$GEO_ID
  meta.4<-meta.4[order(meta.4$TYPE,decreasing = T),]
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB16_mRNA")
  exprMatrix.5 = read.table(file = "HCCDB16_GSE64041_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.5<-data.table::fread(file = "HCCDB16.sample.txt",header = F)
  meta.5<-t(meta.5)%>%as.data.frame()
  colnames(meta.5)<-meta.5[1,]
  meta.5<-meta.5[-1,]
  rownames(meta.5)<-meta.5$GEO_ID
  meta.5<-meta.5[order(meta.5$TYPE,decreasing = T),]
  
  index=intersect(intersect(rownames(exprMatrix.1),rownames(exprMatrix.2)),rownames(exprMatrix.3))
  index=intersect(intersect(index,rownames(exprMatrix.4)),rownames(exprMatrix.5))
  exprMatrix.1 =exprMatrix.1[index,meta.1$GPL6793_GEO_ID]
  exprMatrix.2=exprMatrix.2[index,meta.2$GEO_ID]
  exprMatrix.3=exprMatrix.3[index,meta.3$GEO_ID]
  exprMatrix.4=exprMatrix.4[index,meta.4$GEO_ID]
  exprMatrix.5=exprMatrix.5[index,meta.5$GEO_ID]
  
  tumor.1<-exprMatrix.1[,1:max(which(meta.1$TYPE=="HCC"))]
  normal.1<-exprMatrix.1[,(max(which(meta.1$TYPE=="HCC"))+1):ncol(exprMatrix.1)]
  healthy.2<-exprMatrix.2[,1:max(which(meta.2$TYPE=="Healthy"))]
  tumor.2<-exprMatrix.2[,(max(which(meta.2$TYPE=="Healthy"))+1):max(which(meta.2$TYPE=="HCC"))]
  normal.2<-exprMatrix.2[,(max(which(meta.2$TYPE=="HCC"))+1):ncol(exprMatrix.2)]
  tumor.3<-exprMatrix.3[,1:max(which(meta.3$TYPE=="HCC"))]
  normal.3<-exprMatrix.3[,(max(which(meta.3$TYPE=="HCC"))+1):ncol(exprMatrix.3)]
  tumor.4<-exprMatrix.4[,1:max(which(meta.4$TYPE=="HCC"))]
  normal.4<-exprMatrix.4[,(max(which(meta.4$TYPE=="HCC"))+1):ncol(exprMatrix.4)]
  healthy.5<-exprMatrix.5[,1:max(which(meta.5$TYPE=="Healthy"))]
  tumor.5<-exprMatrix.5[,(max(which(meta.5$TYPE=="Healthy"))+1):max(which(meta.5$TYPE=="HCC"))]
  normal.5<-exprMatrix.5[,(max(which(meta.5$TYPE=="HCC"))+1):ncol(exprMatrix.5)]
  
  exprMatrix<-cbind(exprMatrix.1,exprMatrix.2,exprMatrix.3,exprMatrix.4,exprMatrix.5)
  tumor<-cbind(tumor.1,tumor.2,tumor.3,tumor.4,tumor.5)
  normal<-cbind(normal.1,normal.2,normal.3,normal.4,normal.5)
  healthy<-cbind(healthy.2,healthy.5)
  data<-cbind(healthy.2,healthy.5,normal,tumor)
  data<-cbind(normal,tumor)
  batch=c(rep(c("GSE1","GSE2","GSE3","GSE4","GSE5"),c(97,283,220,168,60)),
          rep(c("GSE1","GSE2","GSE3","GSE4","GSE5"),c(100,268,225,228,60)))
  
  tissue<-rep(c("Normal","Tumor"),c(828,881))
}
{
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB6_mRNA")
  exprMatrix.3 = read.table(file = "HCCDB6_GSE14520-GPL3921_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.3<-data.table::fread(file = "HCCDB6.sample.txt",header = F)
  meta.3<-t(meta.3)%>%as.data.frame()
  colnames(meta.3)<-meta.3[1,]
  meta.3<-meta.3[-1,]
  rownames(meta.3)<-meta.3$GEO_ID
  meta.3<-meta.3[order(meta.3$TYPE,decreasing = T),]
  max(which(meta.3$TYPE=="HCC"))
  
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB13_mRNA")
  exprMatrix.4 = read.table(file = "HCCDB13_GSE63898_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  meta.4<-data.table::fread(file = "HCCDB13.sample.txt",header = F)
  meta.4<-t(meta.4)%>%as.data.frame()
  colnames(meta.4)<-meta.4[1,]
  meta.4<-meta.4[-1,]
  rownames(meta.4)<-meta.4$GEO_ID
  meta.4<-meta.4[order(meta.4$TYPE,decreasing = T),]
  
  index=intersect(rownames(exprMatrix.3),rownames(exprMatrix.4))
  exprMatrix.3=exprMatrix.3[index,meta.3$GEO_ID]
  exprMatrix.4=exprMatrix.4[index,meta.4$GEO_ID]
  
  tumor.3<-exprMatrix.3[,1:max(which(meta.3$TYPE=="HCC"))]
  normal.3<-exprMatrix.3[,(max(which(meta.3$TYPE=="HCC"))+1):ncol(exprMatrix.3)]
  tumor.4<-exprMatrix.4[,1:max(which(meta.4$TYPE=="HCC"))]
  normal.4<-exprMatrix.4[,(max(which(meta.4$TYPE=="HCC"))+1):ncol(exprMatrix.4)]
  
  exprMatrix<-cbind(exprMatrix.3,exprMatrix.4)
  tumor<-cbind(tumor.3,tumor.4)
  normal<-cbind(normal.3,normal.4)
  data<-cbind(normal,tumor)
  batch=c(rep(c("GSE14520-GPL3921","GSE63898-GPL13667"),c(220,168)),
          rep(c("GSE14520-GPL3921","GSE63898-GPL13667"),c(225,228)))
  
  tissue<-rep(c("Normal","Tumor"),c(388,453))
}
###############merge################
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge")
save(tumor,normal,exprMatrix,data,ex_b_sva,batch,tissue,file = "Affy1.Rdata")
load("Affy1.Rdata")
exprMatrix.2 = ex_b_sva
tissue.2<-rep(c("Normal","Tumor"),c(388,453))
exprMatrix.1 = read.table(file = "illumina_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
index=intersect(rownames(exprMatrix.1),rownames(exprMatrix.2))
tissue.1<-rep(c("Normal","Tumor"),c(57+245,181+439))
exprMatrix.1=exprMatrix.1[index,]
exprMatrix.2=exprMatrix.2[index,]
tissue<-c(tissue.1,tissue.2)

data<-cbind(exprMatrix.1,exprMatrix.2)
batch=rep(c("illumina","Affy"),c(922,841))

rt=cbind(ID=rownames(ex_b_limma),as.data.frame(ex_b_limma))
write.table(rt,file = "GPL570_symbol_tpm.txt",quote = F,row.names = F)
#######################xcell#################################
library(xCell)
xcell=xCellAnalysis(exprMatrix)
xcell=xcell[,ann_col$sample]
write.csv(xcell,file = "merge_xcell.csv",quote = F)
xcell.1=xCellAnalysis(exprMatrix.1)
write.csv(xcell.1,file = "illumina_xcell.csv",quote = F)
xcell.2=xCellAnalysis(exprMatrix.2)
write.csv(xcell.2,file = "Affy_xcell.csv",quote = F)
xcell.3=xCellAnalysis(exprMatrix.3)
write.csv(xcell.3,file = "GPL570_xcell.csv",quote = F)
xcell.4=xCellAnalysis(exprMatrix.4)
write.csv(xcell.4,file = "GPL10558_xcell.csv",quote = F)

###
if(F){
  x1=read.csv(file = "GSE121248_xcell.csv",row.names = 1)
  x2=read.csv(file = "GSE41804_xcell.csv",row.names = 1)
  x3=read.csv(file = "HCCDB8_GSE9843_xcell.csv",row.names = 1)
  x4=read.csv(file = "HCCDB4_GSE36376_xcell.csv",row.names = 1)
  x5=read.csv(file = "HCCDB14_GSE43619_xcell.csv",row.names = 1)
  x6=read.csv(file = "HCCDB17_GSE76427_xcell.csv",row.names = 1)
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB6_mRNA")
  x7=read.csv(file = "HCCDB6_GSE14520-GPL3921_xcell.csv",row.names = 1)
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\HCCDB13_mRNA")
  x8=read.csv(file = "HCCDB13_GSE63898_xcell.csv",row.names = 1)
  setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge")
  xcell=cbind(x1,x2,x3,x4,x5,x6,x7,x8)
  write.csv(xcell,file = "merge_xcell.txt")
  exprMatrix = read.table(file = "merge_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
  
  ann_col<-data.frame(sample=colnames(exprMatrix),tissue=tissue,batch=batch)
  ann_col<-ann_col[order(ann_col$tissue),]
  rownames(ann_col)<-ann_col$sample
  xcell<-xcell[,colnames(exprMatrix)]
}

###
source("D:\\Bioinfrolf\\Bioinfrolf\\ICGA-HCC\\JP\\CIBERSORT-2021.R")
rt=cbind(ID=rownames(exprMatrix),as.data.frame(exprMatrix))
write.table(rt,file = "Tumor_symbol_tpm.txt",quote = F,row.names = F)
results0=CIBERSORT("ref.txt", "Tumor_symbol_tpm.txt", perm=1000,QN=TRUE)
results0=t(results0)
results0=data.table::fread("CIBERSORT-Results.txt")
a=results0$Mixture
results0=results0[,2:23]
rownames(results0)=a

ann_col<-data.frame(sample=a,tissue=tissue,batch=batch)
ann_col<-ann_col[order(ann_col$tissue),]
rownames(ann_col)<-ann_col$sample

xcell=t(results0)
colnames(xcell)=a

###
#install_github('ebecht/MCPcounter',ref='master', subdir='Source')
library(MCPcounter)
probesets=read.table("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/MCPcounter-master/Signatures/probesets.txt",sep="\t",stringsAsFactors=FALSE,colClasses="character")
genes=read.table("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/MCPcounter-master/Signatures/genes.txt",
                 sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
results<- MCPcounter.estimate(exprMatrix,featuresType=c("HUGO_symbols")[1],
                              probesets=probesets,
                              genes=genes
)
write.csv(results,file = "merge_MCPcounter.csv")
xcell=results[,ann_col$sample]
save(results,results0,results1,file = "Tumor_immune.Rdata")
load("Tumor_immune.Rdata")

###
GS=read.csv("Tumor_signature_scores_scaled.csv",row.names = 1)
GS_D=Dscaling(GS)
GS_M=MinMax(GS)

####


rownames(results0)
rownames(results1)
rownames(results)

h<-pheatmap::pheatmap(df,scale = "row",
                   show_colnames =F,
                   show_rownames = T,
                   annotation_col = annotation_col,
                   cluster_row = T, 
                   cluster_col = T,
                   border_color ='white',
                   angle_col=45,#fontsize=10,fontsize_col=10,fontsize_row=6,
                   color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                             "RdYlBu")))(100),
                   main ="macrophages")

jpeg(file="macro_score.jpeg",width=10,height=8,units = "in", res = 600)
h
dev.off()

name<-rownames(ex_z)[h$tree_row$order]

mat_MCP<-scale_rows(results)
rownames(mat_MCP)<-rownames(results)
colnames(mat_MCP)<-colnames(results)
mat<-mat[name,]

mat<-as.data.frame(mat)
g_anno2 <- data.frame(SPP1 = ex_z["SPP1",],sample=colnames(ex_z))
rownames(g_anno2)<-colnames(ex_z)
g_anno2<-g_anno2[resn$gene_name,]

g_anno1 <- data.frame(type = c(rep('Normal High Expression',200),rep('TumorCore High Expression',100),
                               rep('TumorEdge High Expression',100),rep('Normal High Expression',45)))
s_anno <- data.frame(value = c("Normal", "TumorEdge","TumorCore"))
s_anno$value <-factor(s_anno$value,levels = c("Normal", "TumorEdge","TumorCore"))
cirmat$variable<-factor(cirmat$variable,levels = c("Normal", "TumorEdge","TumorCore"))

p1 <- ggplot() +
  geom_bar(data = g_anno2,stat = 'identity',
           aes(x = 1:nrow(resn),y = -0.5,fill = importance),
           width = 1,
           color = NA) +
  scale_fill_gradient(low = "#FFEEEE",high = '#FF6767')+
  new_scale("fill") +
  geom_bar(data = g_anno1,stat = 'identity',
           aes(x = 1:nrow(resn),y = -0.25,fill = type),
           width = 1,
           color = NA) +
  scale_fill_manual(name = 'gene annotation',
                    values = 
                      c('Normal High Expression'="#639452",
                        'TumorEdge High Expression'="#F97500",
                        'TumorCore High Expression'="#FC3400"))+
  new_scale("fill") +
  geom_bar(data = s_anno,stat = 'identity',
           aes(x = -2,y = 1,fill = value),
           width = 8,
           color = NA) +
  scale_fill_manual(name = 'Tissue',
                    values = c("Normal"="#00A08A", 
                               "TumorEdge"="#F98400",
                               "TumorCore"="#FF0000"))

p3 <- p1 + 
  new_scale("fill") +
  geom_tile(data = cirmat[which(cirmat$variable == 'Normal'),],
            aes(x = 1:nrow(resn),y = 2.5,fill = value),
            color = 'white') +
  geom_tile(data = cirmat[which(cirmat$variable == 'TumorEdge'),],
            aes(x = 1:nrow(resn),y = 1.5,fill = value),
            color = 'white') +
  geom_tile(data = cirmat[which(cirmat$variable == 'TumorCore'),],
            aes(x = 1:nrow(resn),y = 0.5,fill = value),
            color = 'white') +
  scale_fill_gradient2(midpoint = 0,
                       low = "#5BBCD6",
                       mid = "white",
                       high = "#B40F20") +
  ylim(-3,5)


Genes<-read.table("impMatrix_201.txt",header = T)
index=intersect(Genes$Feature,rownames(exprMatrix.4))
which(index%in%c("MIF","SPP1","CD74","CXCL12"))
Genes$Feature[which(!Genes$Feature%in%index)]


#######################cancer cell####################
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge")
exprMatrix = read.table(file = "seq_all_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
ex<-exprMatrix#log2(as.matrix(exprMatrix)+1)
ex[1:5,1:5]
exprMatrix[1:5,1:5]
ex<-t(ex)
dim(ex)
write.csv(ex,file = "seq_tumor_symbol_tpm.csv")


