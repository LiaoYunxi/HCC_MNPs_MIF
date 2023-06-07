################prepare################
library(tidyverse)
library(stringr)
rm(list = ls())
col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6",
        "#E2D200","#B40F20","#273046","#FD6467")
cols<-c("#5BBCD6","#E2D200","#B40F20","#273046","#FD6467")
col2<-c("#00A08A","#FF0000")
col4<-c("D"="#F98400","F"="#E2D200","IE"="#00A08A","IE/F"="#5BBCD6")
colsurv=c("CT"="#FF0000","AN"="#00A08A","PT"="#F98400")
col5<-c("I"="#D3D3D3","II"="#C89197","III"="#BE505B",
        "IV"="#B40F20","unknown"="white")
col5.1<-c("1"="#D3D3D3","2"="#C89197","3"="#BE505B",
          "4"="#B40F20","unknow"="white")
barplot(1:8, col= col8)

Zscore=function(xcell){
  df<-t(xcell)%>%as.data.frame()
  df<-scale(df,)%>%t()
  return(df)
}

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

Dscaling=function(xcell){
  df<-t(xcell)%>%as.data.frame()
  i1=ceiling(log(max(abs(df[,1])),10))#小数定标的指数
  c1=df[,1]/10^i1
  dat=c1
  for (i in colnames(df)[-1]){
    i1=ceiling(log(max(abs(df[,i])),10))#小数定标的指数
    c1=df[,i]/10^i1
    dat<-cbind(dat,as.data.frame(c1))
  }
  colnames(dat)<-colnames(df)
  rownames(dat)<-rownames(df)
  dat<-t(dat)
  return(dat)
}

MinMax=function(xcell){
  df<-t(xcell)%>%as.data.frame()
  dat=(df[,1]-min(df[,1]))/(max(df[,1])-min(df[,1]))
  for (i in colnames(df)[-1]){
    a=(df[,i]-min(df[,i]))/(max(df[,i])-min(df[,i]))
    dat<-cbind(dat,as.data.frame(a))
  }
  colnames(dat)<-colnames(df)
  rownames(dat)<-rownames(df)
  dat<-t(dat)
  return(dat)
}

CIBERSORT_NAME<-function(x){
  a=x$Mixture
  x=t(x[,2:23])
  colnames(x)=a
  colnames(x)<-gsub("-",".",colnames(x))
  x<-as.data.frame(x)
  return(x)
}

Genes<-read.table("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/table/impMatrix.txt",header = T)
rownames(Genes)<-Genes$Feature

PCAnc<-function(e,name,set,group_list,pcdim){
  dir.create(name)
  library(factoextra)
  e<-t(e)
  e<-as.data.frame(e)
  res.pca <- prcomp(e, scale = TRUE)
  pca<-res.pca$x[,1:pcdim]
  
  jpeg(file=paste0(name,"\\log.feature.pca",set,".jpeg"),width =7,height =5,units = "in", res = 1000)
  fviz_pca_ind(res.pca,
               col.ind = group_list,
               palette =col4,
               pch=15,
               geom = c("point"),
               ellipse.type = "confidence",
               legend.title = "Group",
               repel = TRUE
  )
  dev.off()
  
  data<-data.frame(x=colnames(res.pca$x),y=res.pca$sdev)
  data$x<-factor(data$x,levels = data$x)
  
  jpeg(file=paste0(name,"\\log.pca",pcdim,".elbow.",set,".jpeg"),width =6,height =3,units = "in", res = 1000)
  ggplot(data[1:pcdim,],aes(x=x,y=y))+geom_point(size=2)+
    theme_bw()+theme(panel.grid.major =element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),
                     axis.text.x =  element_text(hjust = 1,vjust = 1,size = 3)
    )+labs(title = sum(data[1:pcdim,2]))
  dev.off()
  
  
  set.seed(123)
  wssplot <- function(data, nc=15, seed=1234){
    wss <- (nrow(data)-1)*sum(apply(data, 2, var))
    for(i in 2:nc){
      set.seed(seed)
      wss[i] <- sum(kmeans(data, centers = i)$withinss)
    }
    plot(1:nc, wss, type = "b",xlab = "Number of Clusters",
         ylab = "Within groups sum of squares")
  }
  
  jpeg(file=paste0(name,"\\log.kmeans.wssplot.",set,".jpeg"),width =6,height =4,units = "in", res = 1000)
  wssplot(pca)
  dev.off()
  
  library(NbClust)
  nb <- NbClust(pca, min.nc = 2, max.nc = 8, method = "kmeans")
  dev.off()
  
  jpeg(file=paste0(name,"\\log.kmeans.NbClust.",set,".jpeg"),width =6,height =4,units = "in", res = 1000)
  barplot(table(nb$Best.n[1,]),
          xlab = "Number of Clusters", ylab = "Number of Criteria",
          main = "Number of Clusters Chosen by 26 Criteria")
  dev.off()
  
  nb<-which.max(as.character(table(nb$Best.n[1,])))
  
  return(nb)
}

Dividgroup<-function(e,name,set,group_list,pcdim,nc,seed){
  cluster3<-CLUSTER3(e,name)
  cluster2<-CLUSTER2(e,name)
  cluster3.w<-CLUSTER3.wt(e)
  cluster2.w<-CLUSTER2.wt(e)
  
  library(factoextra)
  e<-t(e)
  e<-as.data.frame(e)
  res.pca <- prcomp(e, scale = TRUE)
  pca<-res.pca$x[,1:pcdim]
  # 
  # library(mclust)
  # m_clust <- Mclust(as.matrix(pca), G=2:5)
  # summary(m_clust)
  # #str(m_clust)
  # 
  # library(fpc)
  # pamk.best <- pamk(pca)
  # pamk.best$nc
  # library(cluster)
  # pa<-pam(pca,k=nc)
  # 
  # set.seed(seed)
  # kmeans<- kmeans(pca, centers=nc, iter.max = 1000, nstart = 25)
  # # fviz_cluster(kmeans, data = pca)
  # # fviz_cluster(pa, data = pca)
  # #dev.off()
  # 
  library(umap)
  ump <- umap(pca, n_neighbors = 5, random_state = 3)
  library(fpc)
  db <- fpc::dbscan(ump$layout, eps = 0.15, MinPts = 5)
  library(dbscan)
  op<-optics(ump$layout)
  dfump<-as.data.frame(ump$layout)
  dfump$db<-db$cluster
  dfump$op<-op$order
  
  library(Rtsne)
  set.seed(2)
  tsne_out = Rtsne(pca,
                   dims = 2, pca = F,
                   max_iter = 1000,theta = 0.7,
                   perplexity = 20, verbose = F,check_duplicates = F)
  dftsne = as.data.frame(tsne_out$Y)
  library(fpc)
  db <- fpc::dbscan(dftsne, eps = 0.15, MinPts = 5)
  library(dbscan)
  op<-optics(dftsne)
  dftsne$db<-db$cluster
  dftsne$op<-op$order
  
  df<-cbind(dfump,dftsne)
  # 
  # df$fc<-group_list
  # df$km<-kmeans$cluster
  # df$mc<-m_clust$classification
  # df$pa<-pa$clustering
  df$cluster3<-cluster3
  df$cluster2<-cluster2
  df$c3<-cluster3.w
  df$c2<-cluster2.w
  
  return(df)
}

CLUSTER3<-function(e,name){
  e<-Zscore(e)
  d0<-e[rownames(e)%in%G0,]
  d1<-e[rownames(e)%in%G1,]
  d2<-e[rownames(e)%in%G2,]
  
  s0<-apply(d0,2,mean)
  s1<-apply(d1,2,mean)
  s2<-apply(d2,2,mean)
  data<-data.frame(row.names =colnames(d0),TumorEdge=s1,TumorCore=s2
                   ,Normal=s0#,gene=rownames(d0)
  )
  # d0<-colnames(e)[apply(data, 1, which.max)==1]
  # d1<-colnames(e)[apply(data, 1, which.max)==2]
  # d2<-colnames(e)[apply(data, 1, which.max)==3]
  
  cluster<-apply(data, 1, which.max)
  cluster<-ifelse(cluster==1,"AN",ifelse(cluster==2,"PT","CT"))%>%as.factor()
  levels(cluster)<-c("AN","PT","CT")
  return(cluster)
}

CLUSTER2<-function(e,name){
  e<-Zscore(e)
  d0<-e[rownames(e)%in%G0,]
  d1<-e[rownames(e)%in%G1,]
  d2<-e[rownames(e)%in%G2,]
  
  s0<-apply(d0,2,mean)
  s1<-apply(d1,2,mean)
  s2<-apply(d2,2,mean)
  data<-data.frame(row.names =colnames(d0),TumorEdge=s1,TumorCore=s2
  )
  # d0<-colnames(e)[apply(data, 1, which.max)==1]
  # d1<-colnames(e)[apply(data, 1, which.max)==2]
  # d2<-colnames(e)[apply(data, 1, which.max)==3]
  
  cluster<-apply(data, 1, which.max)
  cluster<-ifelse(cluster==2,"PT","CT")%>%as.factor()
  levels(cluster)<-c("AN","PT","CT")
  return(cluster)
}

CLUSTER3.wt<-function(e){
  e<-Zscore(e)
  d0<-e[rownames(e)%in%G0,]
  d1<-e[rownames(e)%in%G1,]
  d2<-e[rownames(e)%in%G2,]
  
  s0<-apply(d0,2,function(x){
    weighted.mean(x,Genes[rownames(d0),2])
  })
  s1<-apply(d1,2,function(x){
    weighted.mean(x,Genes[rownames(d1),2])
  })
  s2<-apply(d2,2,function(x){
    weighted.mean(x,Genes[rownames(d2),2])
  })
  data<-data.frame(row.names =colnames(d0),TumorEdge=s1,TumorCore=s2
                   ,Normal=s0#,gene=rownames(d0)
  )
  # d0<-colnames(e)[apply(data, 1, which.max)==1]
  # d1<-colnames(e)[apply(data, 1, which.max)==2]
  # d2<-colnames(e)[apply(data, 1, which.max)==3]
  
  cluster<-apply(data, 1, which.max)
  cluster<-ifelse(cluster==1,"AN",ifelse(cluster==2,"PT","CT"))%>%as.factor()
  levels(cluster)<-c("AN","PT","CT")
  return(cluster)
}

CLUSTER2.wt<-function(e){
  e<-Zscore(e)
  d0<-e[rownames(e)%in%G0,]
  d1<-e[rownames(e)%in%G1,]
  d2<-e[rownames(e)%in%G2,]
  
  
  s0<-apply(d0,2,function(x){
    weighted.mean(x,Genes[rownames(d0),2])
  })
  s1<-apply(d1,2,function(x){
    weighted.mean(x,Genes[rownames(d1),2])
  })
  s2<-apply(d2,2,function(x){
    weighted.mean(x,Genes[rownames(d2),2])
  })
  
  data<-data.frame(row.names =colnames(d0),TumorEdge=s1,TumorCore=s2
  )
  # d0<-colnames(e)[apply(data, 1, which.max)==1]
  # d1<-colnames(e)[apply(data, 1, which.max)==2]
  # d2<-colnames(e)[apply(data, 1, which.max)==3]
  
  cluster<-apply(data, 1, which.max)
  cluster<-ifelse(cluster==2,"PT","CT")%>%as.factor()
  levels(cluster)<-c("AN","PT","CT")
  return(cluster)
}

UMAP<-function(dfump,set,name,cl,clname,col){
  ggplot(dfump,aes(x=V1,y=V2))+geom_point(aes(color = cl),size=2)+
    theme_bw()+theme(panel.grid.major =element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     title = element_text(size = 15,face = 'bold'),
                     legend.text =element_text(size=10),  # Font size of legend labels.
                     legend.title = element_blank(), 
                     legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )+xlab('UMAP1') + ylab('UMAP2')+
    scale_color_manual(values = col)
  
}

TSNE<-function(dftsne,set,name,cl,clname,col){
  ggplot(dftsne,aes(x=V1,y=V2))+geom_point(aes(color = cl),size=2)+
    theme_bw()+theme(panel.grid.major =element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     title = element_text(size = 15,face = 'bold'),
                     legend.text =element_text(size=10),  # Font size of legend labels.
                     legend.title = element_blank(), 
                     legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )+xlab('tSNE1') + ylab('tSNE2')+scale_color_manual(values = col)
  
}


load("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/scGene.Rdata")
Genes$Groups<-ifelse(Genes$Feature%in%G0,"AN",ifelse(Genes$Feature%in%G1,"PT","CT"))
################all##########
library(Scissor)
load("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/03.MNPs.Rdata")
load("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/ICGC_count.Rdata")
sample<-meta$submitted_sample_id
tissue<-unlist(lapply(sample,function(x){
  unlist(strsplit(x,"_"))[2]
}))
tissue<-ifelse(tissue=="Liver","Normal","Tumor")
table(tissue)
meta$tissue<-tissue
rm(data.raw)
tissue=c(rep(c("Normal","Tumor"),c(50,374)),meta$tissue)
batch=rep(c("TCGA","ICGA"),c(424,445))

load("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/seq.original-1.Rdata")
setwd("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/")
all<-ex_b_sva1
rm(ex_b_sva1)
all[1:5,1:5]

batch.t<-batch[tissue=="Tumor"]
ex.all.t<-all[,tissue=="Tumor"]
ex.ICGC.t<-ex.all.t[,batch.t=="ICGA"]
ex.TCGA.t<-ex.all.t[,batch.t=="TCGA"]

################meta####################################
a=unlist(lapply(colnames(ex.TCGA.t), function(x){
  unlist(paste(unlist(str_split(x,"\\."))[1:3],collapse = "-"))
}))

meta.2<-data.table::fread(file = "./table/HCCDB15_time.txt",header = T)
rownames(meta.2)<-meta.2$Id
meta.2=data.frame(row.names = meta.2$Id,futime=meta.2$futime,fustat=meta.2$fustat)

d<-meta.2[a[which(duplicated(a))],]
rownames(d)<-paste0(a[which(duplicated(a))],c(1:3))
meta.2<-rbind(meta.2,d)
a[which(duplicated(a))]<-paste0(a[which(duplicated(a))],c(1:3))
setdiff(rownames(meta.2),a)
meta.2=meta.2[a,]
colnames(ex.TCGA.t)=a

meta.t=meta[meta$tissue=="Tumor",]%>%as.data.frame()
rownames(meta.t)<-meta.t$icgc_specimen_id
meta.t<-meta.t[colnames(ex.ICGC.t),]

meta.3<-data.table::fread(file = "./table/HCCDB18_time.txt",header = T)
rownames(meta.3)<-meta.3$id
meta.3=data.frame(row.names = meta.3$id,futime=meta.3$futime,fustat=meta.3$fustat)
d<-meta.3[meta.t$icgc_donor_id,]
meta.t<-cbind(meta.t,d)
ex.ICGC.t<-ex.ICGC.t[,meta.t$icgc_specimen_id]

# index<-intersect(Genes$Feature,rownames(ex.ICGC.t))
# ex.ICGC.t<-ex.ICGC.t[index,]
# index<-intersect(Genes$Feature,rownames(ex.TCGA.t))
# ex.TCGA.t<-ex.TCGA.t[index,]

################Scissor####################
batch.t<-batch[tissue=="Tumor"]
head(meta.t)
head(meta.2)
meta.icgc<-meta.t[,8:9]
ex.TCGA.t[1:5,1:5]
ex.ICGC.t[1:5,1:5]

save(meta.t,meta.2,meta.icgc,ex.TCGA.t,ex.ICGC.t,file = "Scissor_meta.Rdata")
names(sce)
bulk_dataset=ex.TCGA.t
bulk_survival=meta.2
bulk_survival$TCGA_patient_barcode=rownames(bulk_survival)
all(colnames(bulk_dataset) == bulk_survival$TCGA_patient_barcode)
paid=intersect(colnames(bulk_dataset),bulk_survival$TCGA_patient_barcode)
bulk_dataset<-bulk_dataset[,paid]
bulk_survival<-bulk_survival[paid,]
phenotype <- bulk_survival[,1:2]
colnames(phenotype) <- c("time", "status")
head(phenotype)

location <- "https://xialab.s3-us-west-2.amazonaws.com/Duanchen/Scissor_data/"
load(url(paste0(location, 'TCGA_LUAD_TP53_mutation.RData')))
head(sce@meta.data)
DimPlot(sce, reduction = 'umap', label = T, label.size = 5)
bulk_dataset[1:5,1:5]

pcadim<-25
set.seed(1)
sce <- sce %>% 
  RunUMAP(reduction = "harmony", dims = 1:pcadim) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pcadim)

sc_dataset<-sce
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.05,
                  family = "cox", Save_file = 'Scissor_HCC_survival.RData')

bulk_dataset=ex.ICGC.t
bulk_survival=meta.icgc
bulk_survival$ICGC_patient_barcode=rownames(bulk_survival)
all(colnames(bulk_dataset) == bulk_survival$TCGA_patient_barcode)
phenotype <- bulk_survival[,1:2]
colnames(phenotype) <- c("time", "status")
head(phenotype)
names(sce)
save(bulk_dataset, sc_dataset, phenotype,file ="ICGC_survival_forScissor.RData")
################################
setwd("/Users/zhaolab3/Documents/scRNA_MNP/Rdata/")
all_xcell=read.csv("./table/all.t_xcell.csv", row.names = 1)
all_MCP<-read.csv("./table/all.t_MCP.csv", row.names = 1)
all_CIBERSORT=data.table::fread("./table/all.t-CIBERSORT-Results.txt")
all_CIBERSORT<-CIBERSORT_NAME(all_CIBERSORT)

seq.mfp.t<-read.csv("./table/seq_log_signature_scores-1.csv",row.names = 1)
seq.mfp.t<-seq.mfp.t[colnames(ex.all.t),]
seq.fc.4<-read.csv("./table/seq-log_tumor_clusters_4-1.csv",row.names = 1)
seq.fc.4<-seq.fc.4[colnames(ex.all.t),]

colnames(seq.mfp.t)#"CAF"
rownames(all_xcell)#"Fibroblasts"
rownames(all_CIBERSORT)
rownames(all_MCP)#"Fibroblasts"
fiber<-data.frame(row.names = colnames(ex.all.t),
                  MFP_CAF=as.numeric(seq.mfp.t[colnames(ex.all.t),"CAF"]),
                  xCell_Fibroblasts=as.numeric(all_xcell["Fibroblasts",colnames(ex.all.t)]),
                  MCP_Fibroblasts=as.numeric(all_MCP["Fibroblasts",colnames(ex.all.t)])
)

ann_col_t<-data.frame(row.names = colnames(ex.all.t),sample=colnames(ex.all.t),batch=batch.t)
ann_col_t$fc4<-seq.fc.4
# rownames(fiber)[1:374]<-a
# ex.all.t<-cbind(ex.TCGA.t,ex.ICGC.t)

ex<-all[rownames(ex.all.t),]
ex["MIF",]->mif
log2(ex.all.t["MIF",]+1)->mif
ex.all.t["MIF",]->mif
head(ann_col_t)
head(fiber)
fiber$MIF=mif
data=na.omit(fiber)
#fiber$MFP_CAF=log(fiber$MFP_CAF+1)
x=5
OPT<-function(x,data){
  ndata<-c(0,0,0,0)%>%as.data.frame()
  data<-data[order(data[,4],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    caf<-mean(data[i:(i+x),1])
    xcell<-mean(data[i:(i+x),2])
    mcp<-mean(data[i:(i+x),3])
    m<-mean(data[i:(i+x),4])
    col<-c(caf,xcell,mcp,m)%>%as.data.frame()
    ndata<-cbind(ndata,col)
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  cdt<-rep(0,6)%>%as.data.frame()
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i,method ="spearman")
    cc<-data.frame(c(cmat[1:3,4],cmat[1:2,3],cmat[2,1]))
    colnames(cc)<-i
    cdt<-cbind(cdt,cc)
    #alist=c(alist,list(ndt.i))
  }
  cdt<-cdt[,2:ncol(cdt)]%>%t()%>%as.data.frame()
  return(cdt)
}
arrydata<-OPTy(50,fiber)
arrydata.op<-OPT(19,fiber)
head(fiber)
arrycdt<-cor(arrydata.op,method = "spearman")
quantile(arrydata.op$MFP_CAF)
quantile(log(arrydata.op$MFP_CAF+1))
quantile(fiber$MFP_CAF)
quantile(fiber$MCP_Fibroblasts)
na.omit(arrycdt)
if(correlation){
  #devtools::install_github ("wilkox/gglmannotate")
  library(gglmannotate)
  jpeg(file="MIF-MCP_Fibroblasts.jpeg",width=4,height=4,units = "in", res = 1000)
  ggplot(data=arrydata.op, aes(x=MIF, y=MCP_Fibroblasts)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (TPM)") + ylab("MCP:Fibroblasts") + 
    scale_colour_manual(values = c("grey"))+
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black', size = 0.5), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.text.y= element_text(size = 10),
          axis.text = element_text(size = 10, color = 'black'), 
          axis.title = element_text(size = 15, color = 'black'),
          axis.text.x = element_text(hjust = 1,vjust = 1,size = 12))+theme(
            panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_smooth(aes(x=MIF, y=MCP_Fibroblasts),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                             place = "topleft")+
    annotate("text",x=100,y=9,label="R:0.5035934",parse=T,size=4)
  dev.off()
  
  jpeg(file="MIF-MFP_CAF.jpeg",width=4,height=4,units = "in", res = 1000)
  ggplot(data=arrydata.op, aes(x=MIF, y=MFP_CAF)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (TPM)") + ylab("MFP_CAF") + 
    scale_colour_manual(values = c("grey"))+
    theme(panel.grid = element_blank(), 
          axis.line = element_line(colour = 'black', size = 0.5), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.text.y= element_text(size = 10),
          axis.text = element_text(size = 10, color = 'black'), 
          axis.title = element_text(size = 15, color = 'black'),
          axis.text.x = element_text(hjust = 1,vjust = 1,size = 12))+theme(
            panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_smooth(aes(x=MIF, y=MCP_Fibroblasts),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                             place = "topleft")+
    annotate("text",x=3.5,y=9,label="R:0.2680815",parse=T,size=4)+ylim(2000,5000)
  dev.off()
  bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
  
  jpeg(file="MIF_CAF.infi.jpeg",width=6,height=6,units = "in", res = 1000)
  pheatmap::pheatmap(na.omit(arrycdt),#scale = "row",
                     #kmeans_k = 3,
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = F, cluster_col = F,
                     border_color ='white',
                     angle_col=45,
                     number_color="black",
                     fontsize_number = 16,
                     display_numbers =round(arrycdt, 2),
                     fontsize=12,
                     #fontsize_col=10,fontsize_row=5,
                     color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                               colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                     legend_breaks=seq(-1,1,0.2),
                     breaks=bk,
                     main ="Fibroblasts") 
  dev.off()
  
}