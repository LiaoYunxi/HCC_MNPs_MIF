library(plyr)
cal_Macro_sig_score <- function(h5ad,Clusters,genes){
  ## h5ad: the result of function parse_h5ad
  ## Clusters : a vector of Clusters
  ## genes : a vector of genes
  
  library(scales)
  genes_filt <- genes[genes %in% rownames(h5ad@assays$RNA@counts)]
  expression_matrix <- h5ad@assays$RNA@counts[genes_filt,rownames(h5ad@meta.data[h5ad@meta.data$ann %in% Clusters,])]
  scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
  colnames(scores) <- "value"
  scores$cluster <- h5ad@meta.data[rownames(scores),]$ann
  average_score <- ddply(scores, .(cluster),function(x){mean(x$value)})
  colnames(average_score) <- c("cluster","sig_score")
  average_score$sig_score_scaled <- rescale(average_score$sig_score,to=c(0,5))
  return(average_score)
}

#### M1 and M2 signatures	##########################
# signature from LZY 
M1.2<-read.csv("D:\\Bioinfrolf\\HCC-SC\\myeloid\\plot\\CellMarker-M1.csv")
M2.2<-read.csv("D:\\Bioinfrolf\\HCC-SC\\myeloid\\plot\\CellMarker-M2.csv")

M1.1 <- c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7')
M2.1 <- c('IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC','VEGFD','EGF','CTSA',
          'CTSB','CTSC','CTSD','TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4')
M1<-c(c(M1.2$Cell.Marker),M1.1)%>%unique()
M2<-c(c(M2.2$Cell.Marker),M2.1,"ARG1","RETN")%>%unique()

a<-c(intersect(M1,M2))
M1<-M1[-which(M1%in%a)]
M2<-M2[-which(M2%in%a)]

{
  C("HLA-DPA1,HLA-DP1A,HLASB,
HLA-DOA,HLA-DNA,HLA-DZA,
HLA-DQB2,HLA-DXB,	
HLA-DQA2,HLA-DXA,
HLA-DPB1,HLA-DP1B,HLA-DRB5,HLA-DRB4,HLA-DRB3,
HLA-DQB1,HLA-DQA1,
HLA-DRA,HLA-DRA1,HLA-DMB,DMB,RING7,HLA-DOB,HLA-DMB,HLA-DMA,DMA,RING6,")
}
M1.tmp<-c("FCGR3A,CD16A,FCG3,FCGR3,IGFR3,FCGR3B,CD16B,FCG3,FCGR3,IGFR3,
FCGR2B,CD32,FCG2,IGFR2,FCGR2C,CD32,FCG2,IGFR2,FCGR2A,CD32,FCG2,FCGR2A1,IGFR2,
FCGR1A,FCG1,FCGR1,IGFR1,
IL7R,IL15RA,IL17RA,IL17R,IL2RA,
IL1R1,IL1R,IL1RA,IL1RT1,
CLEC7A,BGR,CLECSF12,DECTIN1,UNQ539,PRO1082,
SELE,ELAM1,MERTK,MER,IL23A,SGRF")

M1.tmp<-str_split(M1.tmp,",")%>%unlist()
M1.tmp<-str_split(M1.tmp,'\n')%>%unlist()%>%unique()

M2.tmp<-c("RETN,FIZZ3,HXCP1,RSTN,
ARG1,IL10,CXCR4,MSR1,SCARA1,
CLEC7A,BGR,CLECSF12,DECTIN1,CSF1R,FMS,
CLEC4A,CLECSF6,DCIR,LLIR,HDCGC13P,
FCER1A,FCE1A,
CCL24,MPIF2,SCYA24,
TNFSF10,APO2L,TRAIL,VEGF,VTCN1,B7H4")
M2.tmp<-str_split(M2.tmp,",")%>%unlist()
M2.tmp<-str_split(M2.tmp,'\n')%>%unlist()%>%unique()

M1<-intersect(rownames(sce@assays$RNA@data),c(M1,M1.tmp))
M2<-intersect(rownames(sce@assays$RNA@data),c(M2,M2.tmp))
M2<-M2[which(M2!="CD14")]
M1<-M1[order(M1)]
M2<-Immunostimulator
M2<-M2[order(M2)]
paste0(M2,collapse = "、")

M1_res <- data.frame()
M2_res <- data.frame()

cluster <- as.vector(unique(sce@meta.data$ann))
library(plyr)
M1_score <- cal_Macro_sig_score(sce,cluster, M1)
M2_score <- cal_Macro_sig_score(sce,cluster, M2)
M1_res <- rbind(M1_res, M1_score)
M2_res <- rbind(M2_res, M2_score)




#set different color vectors for each interval
col1 = colorRampPalette(c("#e9e9e9", 'red'))(30) #set the order of greys


res <- data.frame(cluster=M1_res$cluster, M1=M1_res$sig_score_scaled, M2=M2_res$sig_score_scaled)
#save(res,file = "scoreM1M2.Rdata")

res$cluster->a
res<-as.matrix(res[,-1])
rownames(res)<-a
# res<-res[c("Mono-VCAN","Mono-FCGR3A","Macr-APOC1","Macr-CXCL10","Macr-APOE","Macr-SPP1","Macr-VCAM1"),]
# colnames(res)<-c("angiogenesis","phagocytosis")
res<-res[c("Macr-APOE","Macr-VCAM1","Macr-RPs","Macr-SPP1","Macr-STAT1","Mono-VCAN","Mono-CD16"),]

res<-res[c("Mono-CD16","Mono-VCAN","Macr-RPs","Macr-STAT1","Macr-SPP1","Macr-APOE","Macr-VCAM1"),]
res<-t(res)
jpeg(file="Geneset scores-1.jpeg",width =6,height =4,units = "in", res = 1000)
pheatmap::pheatmap(res,scale = "none",
                   cellwidth = 32, cellheight = 30,
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = FALSE, cluster_col = FALSE,
                   #display_numbers=T,
                   #number_format="%.2f",number_color="black",
                   angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                   border_color ='white',
                   na_col = "grey",
                   main ="Geneset scores")
dev.off()

#####immune###############
setwd("D:/Bioinfrolf/Bioinfrolf/GSEA/")
immunomodulator<-read.table("D:/Bioinfrolf/Bioinfrolf/GSEA/immunomodulator.txt",header = F)
colnames(immunomodulator)<-c("protein","type","gene")
unique(immunomodulator$type)

chemokine<-immunomodulator$gene[immunomodulator$type=="chemokine"]
receptor<-immunomodulator$gene[immunomodulator$type=="receptor"]
MHC<-immunomodulator$gene[immunomodulator$type=="MHC"]
Immunoinhibitor<-immunomodulator$gene[immunomodulator$type=="Immunoinhibitor"]
Immunostimulator<-immunomodulator$gene[immunomodulator$type=="Immunostimulator"]

cal_Macro_sig_score <- function(h5ad,Clusters,genes){
  library(scales)
  genes_filt <- genes[genes %in% rownames(h5ad@assays$RNA@counts)]
  expression_matrix <- h5ad@assays$RNA@data[genes_filt,rownames(h5ad@meta.data[h5ad@meta.data$ann %in% Clusters,])]
  scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
  colnames(scores) <- "value"
  scores$cluster <- h5ad@meta.data[rownames(scores),]$ann
  average_score <- ddply(scores, .(cluster),function(x){mean(x$value)})
  colnames(average_score) <- c("cluster","sig_score")
  average_score$sig_score_scaled <- rescale(average_score$sig_score,to=c(0,5))
  return(average_score)
}

chemokine_res <- data.frame()
receptor_res <- data.frame()
MHC_res <- data.frame()
Immunoinhibitor_res <- data.frame()
Immunostimulator_res <- data.frame()

cluster <- as.vector(unique(sce@meta.data$ann))
library(plyr)
chemokine_score <- cal_Macro_sig_score(sce,cluster, chemokine)
receptor_score <- cal_Macro_sig_score(sce,cluster, receptor)
MHC_score <- cal_Macro_sig_score(sce,cluster, MHC)
Immunoinhibitor_score <- cal_Macro_sig_score(sce,cluster, Immunoinhibitor)
Immunostimulator_score <- cal_Macro_sig_score(sce,cluster, Immunostimulator)

chemokine_res <- rbind(chemokine_res, chemokine_score)
receptor_res <- rbind(receptor_res, receptor_score)
MHC_res <- rbind(MHC_res, MHC_score)
Immunoinhibitor_res <- rbind(Immunoinhibitor_res, Immunoinhibitor_score)
Immunostimulator_res <- rbind(Immunostimulator_res, Immunostimulator_score)

col1 = colorRampPalette(c("#e9e9e9", 'red'))(30) #set the order of greys

res <- data.frame(cluster=chemokine_res$cluster, chemokine=chemokine_res$sig_score_scaled, 
                  receptor=receptor_res$sig_score_scaled)
res<-data.frame(cluster=chemokine_res$cluster,
                Immunoinhibitor=Immunoinhibitor_res$sig_score_scaled,
                Immunostimulator=Immunostimulator_res$sig_score_scaled,
                MHC=MHC_res$sig_score_scaled)
res1 <- data.frame(cluster=chemokine_res$cluster, chemokine=chemokine_res$sig_score_scaled, 
                   receptor=receptor_res$sig_score_scaled)
save(res,file = "scoreM1M2.Rdata")

res<-res1
res$cluster->a
res<-as.matrix(res[,-1])
rownames(res)<-a
res<-res[c("Mono-CD16","Mono-VCAN","Macr-RPs","Macr-STAT1","Macr-SPP1","Macr-APOE","Macr-VCAM1"),]
res<-t(res)
jpeg(file="Geneset scores-chemokine-4.jpeg",width =6,height = 4,units = "in", res = 1000)
pheatmap::pheatmap(res,scale = "none",
                   #color =colorRampPalette(c("#FEFEC0","#D73027"))(100),
                   #colorRampPalette(c("#46ACC8","#B40F20"))(100) ,
                   #annotation_col = group,
                   cellwidth = 20, cellheight = 20,
                   show_colnames =T,
                   show_rownames = T,
                   cluster_row = FALSE, cluster_col = FALSE,
                   #display_numbers=T,
                   number_format="%.2f",number_color="black",
                   angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                   border_color ='white',
                   na_col = "grey",
                   main ="Geneset scores")
dev.off()