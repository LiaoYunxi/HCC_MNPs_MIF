rm(list = ls())
library(tidyr)
library(tibble)
library(stringr)
library(dplyr)
setwd("/Users/zhaolab3/Documents/scRNA_MNP/RNAseq/")
treat1 <- read.table("4IPP-1.htseq.count.txt",skip = 404,header = F)
treat2 <- read.table("4IPP-2.htseq.count.txt",skip = 424,header = F)
control1 <- read.table("TAM-1.htseq.count.txt",skip = 441,header = F)
control2 <- read.table("TAM-2.htseq.count.txt",skip = 432,header = F)
KD1 <- read.table("SK-1.htseq.count.txt",skip = 405,header = F)
KD2 <- read.table("SK-2.htseq.count.txt",skip = 360,header = F)
NC1 <- read.table("SC-1.htseq.count.txt",skip = 399,header = F)
NC2 <- read.table("SC-2.htseq.count.txt",skip = 422,header = F)

rna<-data.frame(row.names = treat1$V1,gene_id=treat1$V1,
                TAM_1=control1$V2,TAM_2=control2$V2,
                IPP_1=treat1$V2, IPP_2=treat2$V2, 
                NC_1=NC1$V2,NC_2=NC2$V2,
                KD_1=KD1$V2,KD_2=KD2$V2)
head(rna)

humanGTF <- readr::read_delim("humanGTF", delim = "\t",
                              escape_double = FALSE, col_names = FALSE,
                              trim_ws = TRUE) %>% dplyr::select(X1,X3)
colnames(humanGTF) <- c("symbol","gene_id")
humanGTF$gene_id <- str_split(humanGTF$gene_id,"[.]",simplify = T)[,1] # ICGC的基因名字不包括版本号，这里需要对ENSEMBL进行整理，删除“.”和后面的数字
humanGTF <- unique(humanGTF) %>% dplyr::select(gene_id,symbol) # 去重

humanGTF <- humanGTF[which(humanGTF$gene_id %in% rna$gene_id),]

rna1 <- left_join(humanGTF,rna,by="gene_id")
rna2 <- aggregate(x = rna1[,3:ncol(rna2)],
                  by = list(symbol = rna1$symbol),   #按照相同symbol分组，在组内计算
                  FUN = max) %>%   #原文中是计算最大值（max），也可以计算平均值（mean）或者中位数（median）
  column_to_rownames(var = 'symbol')

lib.size<-colSums(rna2)
lib.size
max(lib.size)
min(lib.size)
head(rna2)

data=rna2[,c(1:4)]
group_list=c(rep("TAM",2),rep("4-IPP",2))
#group_list=c(rep("NC",2),rep("SPP1_KD",2))
data=rna2[,]
if(unpperquartile){
  data=data[apply(data,1,mean)>0,]
  dim(data)
  data=data[apply(data[,c(1:2)],1,mean)>0 & apply(data[,c(3,4)],1,mean)>0,]
  dim(data)
  d <- DGEList(counts=data,group=factor(group_list))
  d <- calcNormFactors(d, method = "upperquartile", p = 0.75)
  d$samples$norm.factors
  head(d)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  de.com <- exactTest(d,pair=c("NC","SPP1_KD"))
  results <- topTags(de.com,n = length(data[,1]))
  #write.table(as.matrix(results$table),file="chd6_dge.xls",row.names = TRUE, sep="\t")#输出差异表达基因
  res<-as.matrix(results$table)
  
  df<-as.data.frame(data)
  df$logFC<-res[rownames(df),"logFC"]
  df$PValue<-res[rownames(df),"PValue"]
  
  nrDEG=df
  logFC_cutof = with(nrDEG,mean(abs(nrDEG$logFC)) + 2*sd(abs(nrDEG$logFC)))
  logFC_cutof=0.3#logFC_cutof=0.5
  nrDEG$result = as.factor(ifelse(nrDEG$PValue < 0.05 & abs(nrDEG$logFC) > logFC_cutof,
                                  ifelse(nrDEG$logFC > logFC_cutof ,'UP','DOWN'),'NOT'))
  table(nrDEG$result)
  nrow(nrDEG[nrDEG$result =='UP',]) 
  nrow(nrDEG[nrDEG$result =='DOWN',])
  nrDEG$Label = c(rep(NA,nrow(nrDEG)))
  
  nrDEG<-  nrDEG[order(nrDEG$PValue),]
  up.genes<- head(rownames(nrDEG)[which(nrDEG$result =="UP")],10)
  down.genes<- head(rownames(nrDEG)[which(nrDEG$result =="DOWN")],10)
  
  deg.top10.genes<- c(as.character(up.genes),as.character(down.genes))
  nrDEG$Label[match(deg.top10.genes,rownames(nrDEG))]<- deg.top10.genes
  fdrFilter<-0.05
  
  label<-c("MIF","SPP1","CD74","CD44","STAT3","STAT1","APOC1","APOE")
  nrDEG[label,]
  nrDEG$Label[rownames(nrDEG)=="SPP1"]="SPP1"
  nrDEG$result[rownames(nrDEG)=="SPP1"]="DOWN"
  
  library(ggrepel)
  p2<-ggplot(data=nrDEG, aes(x=logFC, y=-log10(PValue), color=result)) +
    geom_point(alpha=0.4, size=2) +
    theme_set(theme_set(theme_bw(base_size=12)))+
    xlab("log2 fold change") + ylab("-log10 p-value")+
    #ggtitle( this_tile2 ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c("#08519C", "grey", "#A50F15"),limits = c("DOWN","NOT","UP"))+
    theme_bw()+theme(
      panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p2<-p2+ geom_hline(yintercept=-log10(fdrFilter),linetype="dashed",size=0.5)+ 
    geom_vline(xintercept= c(-logFC_cutof,logFC_cutof),linetype="dashed",size=0.5)+
    theme(axis.text.x=element_text(hjust =0.5, vjust =0.5,size=12),
          axis.text.y=element_text(size=10),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12))
  p3<-p2+geom_label_repel(aes(label=nrDEG$Label), 
                          #fontface="bold",
                          box.padding=unit(0.2, "lines"), 
                          point.padding=unit(0.3, "lines"), 
                          # segment.colour = "grey50",
                          na.rm = TRUE)
  pdf("/Users/zhaolab3/Documents/scRNA_MNP/RNAseq/vol.pdf",width =8,height = 5)
  p3
  dev.off()
  
  head(nrDEG)
  write.csv(nrDEG,file = "DEG.csv",quote = F)
  
}

#######GSEA####
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(ggnewscale)
library(tibble)
msigdbr_species() 
hs_msigdbr <- msigdbr(species="Homo sapiens")
hsGO <- msigdbr(species="Homo sapiens",category="C5")
head(hsGO)%>%as.data.frame()

mdeg<-read.csv("/Users/zhaolab3/Documents/scRNA_MNP/RNAseq/MIF_inhibit_DEG.csv",row.names = 1)
logFC_cutof=0.3
mdeg$result = as.factor(ifelse(mdeg$PValue < 0.05 & abs(mdeg$logFC) > logFC_cutof,
                               ifelse(mdeg$logFC > logFC_cutof ,'UP','DOWN'),'NOT'))
deg<-read.csv(paste0(
  "/Users/zhaolab3/Documents/scRNA_MNP/RNAseq/SPP1_KD","_DEG.csv"),row.names = 1)
logFC_cutof=1
deg$result = as.factor(ifelse(deg$PValue < 0.05 & abs(deg$logFC) > logFC_cutof,
                              ifelse(deg$logFC > logFC_cutof ,'UP','DOWN'),'NOT'))

s_dg<-intersect(rownames(deg[deg$result=="DOWN",]),rownames(quan_data))
s_ug<-intersect(rownames(deg[deg$result=="UP",]),rownames(quan_data))
m_dg<-intersect(rownames(mdeg[mdeg$result=="DOWN",]),rownames(quan_data))
m_ug<-intersect(rownames(mdeg[mdeg$result=="UP",]),rownames(quan_data))
sg<-c(s_dg,s_ug)
mg<-c(m_dg,m_ug)

gmt=rbind(data.frame(gs_name=rep("MIF activated genes",length(m_dg)),gene_symbols=m_dg),
          data.frame(gs_name=rep("MIF repressed genes",length(m_ug)),gene_symbols=m_ug),
          data.frame(gs_name=rep("SPP1 activated genes",length(s_dg)),gene_symbols=s_dg),
          data.frame(gs_name=rep("SPP1 repressed genes",length(s_ug)),gene_symbols=s_ug))
gmt=as.tibble(gmt)

set.seed(100)
dt=mdeg
dt$score=dt$logFC
scale(dt$score)
dt=arrange(dt,desc(score))
geneList = dt$score
geneList[which(geneList==Inf)]=1000
which(is.na(geneList))
names(geneList) = as.character(rownames(dt))
geneList=sort(geneList,decreasing = T)
head(geneList)

gseaRes <- GSEA(geneList = geneList,
                TERM2GENE = gmt,
                minGSSize    = 15,
                maxGSSize    = 2000,
                pvalueCutoff = "none",
                pAdjustMethod = "none",
                verbose      = FALSE,
                #scoreType = "pos",
                by="fgsea")

data_ga <- data.frame(gseaRes) %>%mkGSEA()#%>%filter(pvalue < 0.05)

