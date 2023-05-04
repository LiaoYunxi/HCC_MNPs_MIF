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

Genes<-read.table("D:/Bioinfrolf/Bioinfrolf/model/impMatrix.txt",header = T)
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


load("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/scGene.Rdata")
Genes$Groups<-ifelse(Genes$Feature%in%G0,"AN",ifelse(Genes$Feature%in%G1,"PT","CT"))
#write.csv(Genes,file = "stable4.csv",quote = F)
################seq-all#######################
load("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/seq/ICGC_count.Rdata")
sample<-meta$submitted_sample_id
tissue<-unlist(lapply(sample,function(x){
  unlist(strsplit(x,"_"))[2]
}))
tissue<-ifelse(tissue=="Liver","Normal","Tumor")
table(tissue)
meta$tissue<-tissue
rm(data.raw)

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")
tissue=c(rep(c("Normal","Tumor"),c(50,374)),meta$tissue)
batch=rep(c("TCGA","ICGA"),c(424,445))

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")

load("seq.original-1.Rdata")

all<-ex_b_sva1
rm(ex_b_sva1)

all_xcell=read.csv(".\\results\\all.t_xcell.csv", row.names = 1)
all_MCP<-read.csv(".\\results\\all.t_MCP.csv", row.names = 1)
all_CIBERSORT=data.table::fread(".\\results\\all.t-CIBERSORT-Results.txt")
all_CIBERSORT<-CIBERSORT_NAME(all_CIBERSORT)

batch.t<-batch[tissue=="Tumor"]
ex.all.t<-all[,tissue=="Tumor"]
ex.ICGC.t<-ex.all.t[,batch.t=="ICGA"]
ex.TCGA.t<-ex.all.t[,batch.t=="TCGA"]

seq.mfp.t<-read.csv(".\\results\\seq_log_signature_scores-1.csv",row.names = 1)
seq.mfp.t<-seq.mfp.t[colnames(ex.all.t),]
seq.fc.4<-read.csv(".\\results\\seq-log_tumor_clusters_4-1.csv",row.names = 1)
seq.fc.4<-seq.fc.4[colnames(ex.all.t),]

if(Normal+tumor){
  e<-all
  group_list<-tissue#rep(c("Normal","Tumor"),c(50,374))#meta$tissue#
  colnames(e)<-group_list
  
  library(factoextra)
  e<-t(e)
  e<-as.data.frame(e)
  res.pca <- prcomp(e, scale = TRUE)
  jpeg(file="feature_pca_noMOD_tissue.jpeg",width =7,height =5,units = "in", res = 1000)
  fviz_pca_ind(res.pca,
               axes = c(1,2),
               col.ind = group_list,
               palette =col3,
               pch=15,
               geom = c("point"),
               #addEllipses = TRUE, # Concentration ellipses
               ellipse.type = "confidence",
               legend.title = "Group",
               repel = TRUE
  )
  dev.off()
  
  library(scatterplot3d)
  data_3d = res.pca$x[,1:4]
  data_3d <- data.frame(data_3d)
  data_3d[1:5,1:4]
  data_3d$type <- batch
  c(table(data_3d$type))
  mycolour <- c(rep(col3[1:2],c(table(data_3d$type))))
  mycolour<-ifelse(tissue=="Tumor","#FF0000","#00A08A")
  jpeg(file="feature_pca_all_noMOD_tissue.jpeg",width =7,height =5,units = "in", res = 1000)
  scatterplot3d(data_3d$PC1,data_3d$PC2,data_3d$PC3,
                color = mycolour, # 样本点颜色
                pch =16, #c(rep(15,50), rep(16,50), rep(17,50)), # 样本点类型
                cex.symbols = 0.5, # 样本点的大小
                font.lab = 2,   #标签大小
                font.axis = 2,
                xlab = 'PC1(20.8%)', # x轴标签
                ylab = 'PC2 (4.9%)',
                zlab = 'PC3 (4.7%)',
                main = '3D-PCA', # 标题
                mar=c(3,2.5,3,1.5)+0.1)
  dev.off()
  library(umap)
  ump<- umap(t(all), n_neighbors = 15, random_state = 123)
  dfump<-as.data.frame(ump$layout)
  # ump.pca<- umap(res.pca$x, n_neighbors = 15, random_state = 123)
  # dfump<-as.data.frame(ump.pca$layout)
  library(Rtsne)
  set.seed(123) # 设置随机数种子
  ex.all<-all[!duplicated(rownames(all)),!duplicated(colnames(all))]
  ex.all<-all[!duplicated(all),]
  d<-all[duplicated(all),]
  tsne_out = Rtsne(t(ex.all),
                   dims = 2, pca = T, 
                   max_iter = 1000,theta = 0.7,
                   perplexity = 20, verbose = F,check_duplicates=FALSE)
  dftsne = as.data.frame(tsne_out$Y)
  
  jpeg(file="gene_umap_all_noMOD.jpeg",width =6,height =5,units = "in", res = 1000)
  ggplot(dfump,aes(x=V1,y=V2))+geom_point(aes(color = batch),size=2)+
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
    scale_color_manual(values = col3)
  dev.off()
  
  jpeg(file="gene_tsne_all_noMOD_tissue.jpeg",width =6,height =5,units = "in", res = 1000)
  ggplot(dftsne,aes(x=V1,y=V2))+geom_point(aes(color = tissue),size=2)+
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
    )+xlab('tSNE1') + ylab('tSNE2')+
    scale_color_manual(values = col2)
  dev.off()
}

macro<-data.frame(row.names = colnames(ex.all.t),
                  MFP_M1=as.numeric(seq.mfp.t[colnames(ex.all.t),"M1_signatures"]),
                  MFP_M=as.numeric(seq.mfp.t[colnames(ex.all.t),"Macrophages"]),
                  xCell_Monocytes=as.numeric(all_xcell["Monocytes",colnames(ex.all.t)]),
                  xCell_M0=as.numeric(all_xcell["Macrophages",colnames(ex.all.t)]),
                  xCell_M1=as.numeric(all_xcell["Macrophages M1",colnames(ex.all.t)]),
                  xCell_M2=as.numeric(all_xcell["Macrophages M2",colnames(ex.all.t)]),
                  CIBERSORT_Monocytes=as.numeric(all_CIBERSORT["Monocytes",colnames(ex.all.t)]),
                  CIBERSORT_M0=as.numeric(all_CIBERSORT["Macrophages M0",colnames(ex.all.t)]),
                  CIBERSORT_M1=as.numeric(all_CIBERSORT["Macrophages M1",colnames(ex.all.t)]),
                  CIBERSORT_M2=as.numeric(all_CIBERSORT["Macrophages M2",colnames(ex.all.t)]),
                  MCP_Monocytic=as.numeric(all_MCP["Monocytic lineage",colnames(ex.all.t)])
)
cd8t<-data.frame(row.names = colnames(ex.all.t),
                 MFP_T_cell_traffic=as.numeric(seq.mfp.t[colnames(ex.all.t),"T_cell_traffic"]),
                 MFP_T_cells=as.numeric(seq.mfp.t[colnames(ex.all.t),"T_cells"]),
                 xCell_CD8_T=as.numeric(all_xcell["CD8+ T-cells",colnames(ex.all.t)]),
                 xCell_CD8_Tem=as.numeric(all_xcell["CD8+ Tem",colnames(ex.all.t)]),
                 xCell_CD8_Tn=as.numeric(all_xcell["CD8+ naive T-cells",colnames(ex.all.t)]),
                 xCell_CD8_Tcm=as.numeric(all_xcell["CD8+ Tcm",colnames(ex.all.t)]),
                 CIBERSORT_CD8_T=as.numeric(all_CIBERSORT["T cells CD8",colnames(ex.all.t)]),
                 MCP_CD8_T=as.numeric(all_MCP["CD8 T cells",colnames(ex.all.t)]))

a=unlist(lapply(colnames(ex.TCGA.t), function(x){
  unlist(paste(unlist(str_split(x,"\\."))[1:3],collapse = "-"))
}))
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\prognosis")
meta.2<-data.table::fread(file = "HCCDB15_time.txt",header = T)
rownames(meta.2)<-meta.2$Id
meta.2=data.frame(row.names = meta.2$Id,futime=meta.2$futime,fustat=meta.2$fustat)

d<-meta.2[a[which(duplicated(a))],]
rownames(d)<-paste0(a[which(duplicated(a))],c(1:3))
meta.2<-rbind(meta.2,d)
a[which(duplicated(a))]<-paste0(a[which(duplicated(a))],c(1:3))
setdiff(rownames(meta.2),a)
meta.2=meta.2[a,]
colnames(ex.TCGA.t)=a

meta.t=meta[meta$tissue=="Tumor",]
meta.t<-meta.t[meta.t$icgc_specimen_id%in%colnames(ex.ICGC.t),]

meta.3<-data.table::fread(file = "HCCDB18_time.txt",header = T)
rownames(meta.3)<-meta.3$id
meta.3=data.frame(row.names = meta.3$id,futime=meta.3$futime,fustat=meta.3$fustat)
d<-meta.3[meta.t$icgc_donor_id,]
meta.t<-cbind(meta.t,d)

ex.ICGC.t<-ex.ICGC.t[,meta.t$icgc_specimen_id]
ex.all.t<-cbind(ex.TCGA.t,ex.ICGC.t)
dim(all)
dim(ex.all.t)
243+374

ann_col_t<-data.frame(row.names = colnames(ex.all.t),sample=colnames(ex.all.t),batch=batch.t)
ann_col_t$fc4<-seq.fc.4
rownames(macro)[1:374]<-a
rownames(cd8t)[1:374]<-a

index<-intersect(Genes$Feature,rownames(ex.all.t))
ex.all.t<-ex.all.t[index,]

e<-ex.all.t
name<-"all.tumor"
set="seq"
group_list<-ann_col_t$fc4
pcdim<-40

nc<-PCAnc(e,name,set,group_list,pcdim)

nc<-2
seed=1
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
  
  library(mclust)
  m_clust <- Mclust(as.matrix(pca), G=2:5)
  summary(m_clust)
  #str(m_clust)
  
  library(fpc)
  pamk.best <- pamk(pca)
  pamk.best$nc
  library(cluster)
  pa<-pam(pca,k=nc)
  
  set.seed(seed)
  kmeans<- kmeans(pca, centers=nc, iter.max = 1000, nstart = 25)
  # fviz_cluster(kmeans, data = pca)
  # fviz_cluster(pa, data = pca)
  #dev.off()
  
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
  
  df$fc<-group_list
  df$km<-kmeans$cluster
  df$mc<-m_clust$classification
  df$pa<-pa$clustering
  df$cluster3<-cluster3
  df$cluster2<-cluster2
  df$c3<-cluster3.w
  df$c2<-cluster2.w
  
  return(df)
}

df<-Dividgroup(e,name,set,group_list,pcdim,nc,seed)

dfump<-df[,c(1:4,9:ncol(df))]
dftsne<-df[,5:ncol(df)]

cl<-as.factor(dfump[,5])
clname<-colnames(dfump)[5]
jpeg(file=paste0(name,"\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
UMAP(dfump,set,name,cl,clname,col=col4)
dev.off()
jpeg(file=paste0(name,"\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
TSNE(dftsne,set,name,cl,clname,col = col4)
dev.off()

i=6
i=7
i=8

for(i in 6:8){
  cl<-as.factor(dfump[,i])
  clname<-colnames(dfump)[i]
  jpeg(file=paste0(name,"\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  UMAP(dfump,set,name,cl,clname,col=cols)
  dev.off()
  jpeg(file=paste0(name,"\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  TSNE(dftsne,set,name,cl,clname,col = cols)
  dev.off()
}

i=9
i=10
i=11
i=12

for(i in 9:12){
  cl<-as.factor(dfump[,i])
  clname<-colnames(dfump)[i]
  jpeg(file=paste0(name,"\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  UMAP(dfump,set,name,cl,clname,col=colsurv)
  dev.off()
  jpeg(file=paste0(name,"\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  TSNE(dftsne,set,name,cl,clname,col = colsurv)
  dev.off()
}

#load("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/merge-N/intergrated.Rdata")

df<-cbind(dfump,macro,cd8t)
df$rank1=ifelse(df$CIBERSORT_M2>=median(df$CIBERSORT_M2),"High","Low")
df$rank2=ifelse(df$xCell_M2>=median(df$xCell_M2),"High","Low")
df$rank3=ifelse(df$CIBERSORT_M1>=median(df$CIBERSORT_M1),"High","Low")
df$rank4=ifelse(df$xCell_M1>=median(df$xCell_M1),"High","Low")
df$rank5=ifelse(df$MFP_M >=median(df$MFP_M),"High","Low")
df$rank6=ifelse(df$MFP_M1 >=median(df$MFP_M1),"High","Low")
table(df$rank1,df$c3)

if(vn){
  library(UpSetR)   
  library(gplots)
  vn<-data.frame(row.names = rownames(df),
                 CIBERSORT_M2=ifelse(df$rank1=="High",1,0),
                 xCell_M2=ifelse(df$rank2=="High",1,0),
                 CIBERSORT_M1=ifelse(df$rank3=="High",1,0),
                 xCell_M1=ifelse(df$rank4=="High",1,0))
  setcol<-c("CIBERSORT_M2"="#B40F20","xCell_M2"="#B40F20",
            "CIBERSORT_M1"="#5BBCD6","xCell_M1"="#5BBCD6")
  
  upset(vn, nsets = 4, sets = c("CIBERSORT_M2","xCell_M2","CIBERSORT_M1","xCell_M1"),
        mb.ratio = c(0.6,0.4),#调整上下两部分的比例
        order.by = c("degree"),#keep.order=T,
        sets.bar.color = setcol,point.size=4,line.size=1,text.scale = 1.5,
        queries = list(                     list(query = intersects,
                                                 params = list("CIBERSORT_M2","xCell_M2"),
                                                 color ="#B40F20", active = T),
                                            list(query = intersects,
                                                 params = list("CIBERSORT_M1","xCell_M1"), 
                                                 color ="#5BBCD6", active = T)))
}

cl=c()
for(i in 1:ncol(df)){
  a=class(df[,i])
  cl=c(cl,a)
}

a=chisq.test(df$fc,df$c3)
a$p.value
a=c(a$statistic,a$p.value)
chisq=data.frame(a)
for(i in colnames(df)[cl!="numeric"]){
  a=chisq.test(df$fc,df[,i])
  a=c(a$statistic,a$p.value)
  chisq=cbind(chisq,data.frame(a))
}
chisq=chisq[,-1]
colnames(chisq)=colnames(df)[cl!="numeric"]

if(umap){
  jpeg(file="feature_umap_all_CIBERSORT_M1.jpeg",width =6,height =5,units = "in", res = 1000)
  ggplot(df,aes(x=V1,y=V2))+geom_point(aes(color = rank(CIBERSORT_M1)),size=2)+
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
    )+xlab('UMAP1') + ylab('UMAP2')+scale_color_gradient(low = "white", high = "#B40F20")
  dev.off()
}

library(ggalluvial)
data=df[,c("fc","rank3","rank4","rank1","rank2","rank5","rank6","c3","c2")]
data$fc=factor(data$fc,levels = c("D","F","IE/F","IE"))
data$fcad=factor(data$fcad,levels = c("D","F","IE/F","IE"))
data$rank4=factor(data$rank4,levels = c("Low","High"))
data$rank3=factor(data$rank3,levels = c("Low","High"))
data$rank1=factor(data$rank1,levels = c("Low","High"))
data$rank2=factor(data$rank2,levels = c("Low","High"))
levels(data$c3)<-c("CT","PT","AN")
levels(df$cluster3)<-c("CT","PT","AN")
unique(data$rank3)

if(sankey){
  jpeg(file="feature_all_sankey1.jpeg",width=6,height=4,units = "in", res = 600)
  ggplot(data = data,
         aes(axis1=c3,axis2= fc,label = after_stat(stratum))) +
    scale_x_discrete(limits = c("sample group", "MFP"), expand = c(.01, .05)) +
    geom_alluvium(aes(fill = c3)) +
    geom_stratum() + geom_text(stat = "stratum") +
    theme_minimal() +theme(panel.grid = element_blank())+
    scale_fill_manual(values =colsurv)
  dev.off()
  
  jpeg(file="fc_Proportion.jpeg",width=6,height=4,units = "in", res = 2000)
  ggplot(data, aes(factor(c3)))+ 
    geom_bar(aes(fill = fc), position = "fill")+ xlab("")+
    ylab("Proportion")+theme(legend.title=element_blank(),
                             panel.grid.major =element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.background = element_blank(),
                             panel.border = element_blank(),
                             strip.background.x = element_blank())+ 
    scale_fill_manual(values = c("#F98400","#E2D200","#00A08A","#5BBCD6"))+
    theme(axis.text.x=element_text(hjust=1))
  dev.off()
}

ex<-all[rownames(ex.all.t),]
ex["CD74",]->cd74
ex["SPP1",]->spp1
ex["MIF",]->mif

arrydata<-rbind(MIF=mif,SPP1=spp1,CD74=cd74)%>%t()%>%as.data.frame()
arrydata<-log2(arrydata+1)%>%na.omit()%>%as.data.frame()
cor(arrydata[,1:2])

OPT<-function(x,data){
  ndata<-c(0,0,0)%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    s<-mean(data[i:(i+x),1])
    m<-mean(data[i:(i+x),2])
    cd<-mean(data[i:(i+x),3])
    col<-c(s,m,cd)%>%as.data.frame()
    ndata<-cbind(ndata,col)
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  #alist<-list()
  cdt<-c(0,0,0)%>%as.data.frame()
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-data.frame(c(cmat[1,2:3],cmat[2,3]))
    colnames(cc)<-i
    cdt<-cbind(cdt,cc)
    #alist=c(alist,list(ndt.i))
  }
  cdt<-cdt[,-1]%>%t()%>%as.data.frame()
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(4,arrydata)

arrydata<-cbind(data.frame(SPP1=spp1,MIF=mif),macro,cd8t)%>%as.data.frame()

OPT<-function(x,data){
  ndata<-c(rep(0,ncol(data)))%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    col<-c()
    for(j in 1:ncol(data)){
      s<-mean(data[i:(i+x),j])
      col<-c(col,s)
    }
    ndata<-cbind(ndata,as.data.frame(col))
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  cdt<-c(rep(0,ncol(data)*2-4))%>%as.data.frame()%>%t()
  colnames(cdt)<-c(paste("SPP1",colnames(data)[-2:-1],sep = "-"),paste("MIF",colnames(data)[-2:-1],sep = "-"))
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-t(data.frame(c(cmat[1,3:ncol(data)],cmat[2,3:ncol(data)])))
    colnames(cc)<-colnames(cdt)
    cdt<-rbind(cdt,cc)
  }
  cdt<-cdt[-1,]
  rownames(cdt)<-c(2:(y+1))
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(9,arrydata)

ce<-cor(t(ex))

if(correlation){
  library(gglmannotate)
  jpeg(file=paste0(name,"\\seq-",set,"-MIF-SPP1.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data=arrydata.op, aes(x=MIF, y=SPP1)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (logTPM)") + ylab("SPP1 (logTPM)") + 
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
    geom_smooth(aes(x=MIF, y=SPP1),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                             place = "topleft")+
    annotate("text",x=3.5,y=9,label="R:0.68",parse=T,size=4)
  dev.off()
  
  jpeg(file=paste0(name,"\\cor.all.intergrated.infi.jpeg"),width=8,height=6,units = "in", res = 1000)
  pheatmap::pheatmap(na.omit(arrycdt),#scale = "row",
                     #kmeans_k = 3,
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = F, cluster_col = T,
                     border_color ='white',
                     angle_col=45,
                     number_color="gray",
                     display_numbers =round(arrycdt, 2),
                     fontsize=5,
                     #fontsize_col=10,fontsize_row=5,
                     color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                               colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                     legend_breaks=seq(-1,1,0.2),
                     breaks=bk,
                     main ="correlation optimization") 
  dev.off()
}

if(circle){
  cora<-rbind(SPP1=arrycdt[9,1:(ncol(arrycdt)/2)],MIF=arrycdt[9,(ncol(arrycdt)/2+1):ncol(arrycdt)])%>%as.data.frame()
  colnames(cora)<-colnames(arrydata)[c(-2:-1)]
  h<-pheatmap::pheatmap(t(cora),scale = "none",
                        show_colnames =T,
                        show_rownames = T,
                        cluster_row = T, cluster_col = F,
                        border_color =NA,
                        number_color="black",
                        display_numbers =round(t(cora), 2),
                        cellwidth = 30, cellheight = 30,
                        angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                        na_col = "grey",
                        color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                                  colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                        legend_breaks=seq(-1,1,0.5),
                        breaks=bk,
                        main ="correlation") 
  jpeg(file=paste0(name,"\\cor.all.intergrated.jpeg"),width=6,height=8,units = "in", res = 1000)
  h
  dev.off()
  
  fdf<-t(cora)
  name<-rownames(fdf)[h$tree_row$order]
  fdf<-fdf[name,]
  
  library(reshape2)
  cirmat<-melt(fdf)
  head(cirmat)
  
  resn <- cirmat %>% filter(Var2 == 'MIF')
  rese <- cirmat %>% filter(Var2 == 'SPP1')
  
  resn$ang <- seq(from = (360/nrow(resn)) / 1.5,
                  to = (1.5* (360/nrow(resn))) - 360,
                  length.out = nrow(resn)) + 80
  
  resn$hjust <- 0
  resn$hjust[which(resn$ang < -90)] <- 1
  resn$ang[which(resn$ang < -90)] <- (180+resn$ang)[which(resn$ang < -90)]
  rese$hjust<-resn$hjust
  rese$ang<-resn$ang
  cirmat$var <- rep(1:nrow(resn),2)
  
  s_anno <- data.frame(value = c("SPP1","MIF"))
  s_anno$value <-factor(s_anno$value,levels = c("SPP1","MIF"))
  cirmat$Var2<-factor(cirmat$Var2,levels = c("SPP1","MIF"))
  
  library(gplots)
  library(ggdendro)
  library(ggh4x)
  library(ggnewscale)
  library(ggrepel)
  p1 <- ggplot() +
    geom_bar(data = s_anno,stat = 'identity',
             aes(x = 0,y = 1,fill = value),
             width = 1,
             color = NA) +
    scale_fill_manual(name = 'Gene',
                      values = c("MIF"="#E2D200", 
                                 "SPP1"="#F98400"))
  p2<-p1+ new_scale("fill")+
    geom_tile(data = cirmat[which(cirmat$Var2 == 'MIF'),],
              aes(x = 1:nrow(resn),y = 0.5,fill = value),
              color = 'white') +
    geom_tile(data = cirmat[which(cirmat$Var2 == 'SPP1'),],
              aes(x = 1:nrow(resn),y = 1.5,fill = value),
              color = 'white')+  scale_fill_gradient2(midpoint = 0,
                                                      low = "#5BBCD6",
                                                      mid = "white",
                                                      high = "#B40F20") +
    ylim(-1,5)+xlim(-1,20)
  
  jpeg(file=paste0(name,"\\cor.all.intergrated.circle.jpeg"),width=15,height=10,units = "in", res = 1000)
  p2 + coord_polar(theta = 'x') +theme_void() +new_scale("colour")+
    geom_text(data = resn,
              aes(x = as.numeric(rownames(resn)),
                  y = 2.2,
                  label = as.character(Var1), angle = ang, hjust = hjust),
              size = 4)+
    geom_text(data = resn,
              aes(x = as.numeric(rownames(resn)),
                  y = c(0.25),
                  label = round(value,2), angle = ang, hjust = hjust),
              size = 3.5)+
    geom_text(data = rese,
              aes(x = as.numeric(rownames(rese)),
                  y = c(1.25),
                  label = round(value,2),angle = ang, hjust = hjust),
              size = 4.5)#+NoLegend()
  # annotate("text", x = as.numeric(rownames(rese)), 
  #          y = 0.25, label = round(rese$value,2), color = "white", size = 4)
  dev.off()
}

df.infi<-cbind(macro,cd8t)

annotation_col = data.frame(
  #MFP=factor(df$fc),
  #batch=factor(batch.t)
  sg=factor(df$c3)
)

rownames(annotation_col) =colnames(ex.all.t)

jpeg(file=paste0(name,"\\cor.all.intergrated.infi.jpeg"),width=8,height=6,units = "in", res = 1000)
pheatmap::pheatmap(na.omit(macro),scale = "column",
                   #kmeans_k = 3,
                   show_colnames =T,
                   show_rownames = F,
                   cluster_row = T, cluster_col =T,
                   annotation_row = annotation_col,
                   border_color ='white',
                   angle_col=45,
                   fontsize=5,
                   #fontsize_col=10,fontsize_row=5,
                   # color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                   #           colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                   #legend_breaks=seq(-1,1,0.2),
                   #breaks=bk,
                   main ="correlation optimization") 
dev.off()

################seq-divide################
conNum=50           
treatNum=374
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")
exprMatrix.2 = read.table(file = "TCGA_symbol_tpm.txt",
                          header=TRUE,row.names=1, as.is=TRUE)

exprMatrix.3 = read.table(file = "ICGC_symbol_tpm.txt",
                          header=TRUE,row.names=1, as.is=TRUE)
exprMatrix.3[1:5,1:5]

ex.TCGA.t<-exprMatrix.2[,51:ncol(exprMatrix.2)]
ex.ICGC.t<-exprMatrix.3[,meta$tissue=="Tumor"]

ICGC_xcell=read.csv(".\\results\\ICGC_xcell.csv", row.names = 1)
ICGC_MCP<-read.csv(".\\results\\ICGC_MCP.csv", row.names = 1)
ICGC_CIBERSORT=data.table::fread(".\\results\\ICGC-CIBERSORT-Results.txt")
ICGC_CIBERSORT<-CIBERSORT_NAME(ICGC_CIBERSORT)
TCGA_xcell=read.csv(".\\results\\TCGA_xcell.csv", row.names = 1)
TCGA_MCP<-read.csv(".\\results\\TCGA_MCP.csv", row.names = 1)
TCGA_CIBERSORT=data.table::fread(".\\results\\TCGA-CIBERSORT-Results.txt")
TCGA_CIBERSORT<-CIBERSORT_NAME(TCGA_CIBERSORT)

TCGA.mfp.t<-read.csv("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq\\results\\TCGA_log_signature_scores_scaled.csv",row.names = 1)
ICGC.mfp.t<-read.csv("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq\\results\\ICGC_log_signature_scores_scaled.csv",row.names = 1)
TCGA.fc.4<-read.csv(".\\results\\TCGA-log_tumor_clusters_4.csv",row.names = 1)
ICGC.fc.4<-read.csv(".\\results\\ICGC-log_tumor_clusters_4.csv",row.names = 1)

if(score){
  macro.TCGA<-data.frame(row.names = colnames(ex.TCGA.t),
                         MFP_M1=as.numeric(TCGA.mfp.t[colnames(ex.TCGA.t),"M1_signatures"]),
                         MFP_M=as.numeric(TCGA.mfp.t[colnames(ex.TCGA.t),"Macrophages"]),
                         xCell_Monocytes=as.numeric(TCGA_xcell["Monocytes",colnames(ex.TCGA.t)]),
                         xCell_M0=as.numeric(TCGA_xcell["Macrophages",colnames(ex.TCGA.t)]),
                         xCell_M1=as.numeric(TCGA_xcell["Macrophages M1",colnames(ex.TCGA.t)]),
                         xCell_M2=as.numeric(TCGA_xcell["Macrophages M2",colnames(ex.TCGA.t)]),
                         CIBERSORT_Monocytes=as.numeric(TCGA_CIBERSORT["Monocytes",colnames(ex.TCGA.t)]),
                         CIBERSORT_M0=as.numeric(TCGA_CIBERSORT["Macrophages M0",colnames(ex.TCGA.t)]),
                         CIBERSORT_M1=as.numeric(TCGA_CIBERSORT["Macrophages M1",colnames(ex.TCGA.t)]),
                         CIBERSORT_M2=as.numeric(TCGA_CIBERSORT["Macrophages M2",colnames(ex.TCGA.t)]),
                         MCP_Monocytic=as.numeric(TCGA_MCP["Monocytic lineage",colnames(ex.TCGA.t)])
  )
  cd8t.TCGA<-data.frame(row.names = colnames(ex.TCGA.t),
                        MFP_T_cell_traffic=as.numeric(TCGA.mfp.t[colnames(ex.TCGA.t),"T_cell_traffic"]),
                        #MFP_T_cells=as.numeric(TCGA.mfp.t[colnames(ex.TCGA.t),"T_cells"]),
                        xCell_CD8_T=as.numeric(TCGA_xcell["CD8+ T-cells",colnames(ex.TCGA.t)]),
                        xCell_CD8_Tem=as.numeric(TCGA_xcell["CD8+ Tem",colnames(ex.TCGA.t)]),
                        xCell_CD8_Tn=as.numeric(TCGA_xcell["CD8+ naive T-cells",colnames(ex.TCGA.t)]),
                        xCell_CD8_Tcm=as.numeric(TCGA_xcell["CD8+ Tcm",colnames(ex.TCGA.t)]),
                        CIBERSORT_CD8_T=as.numeric(TCGA_CIBERSORT["T cells CD8",colnames(ex.TCGA.t)]),
                        MCP_CD8_T=as.numeric(TCGA_MCP["CD8 T cells",colnames(ex.TCGA.t)]))
}

a=unlist(lapply(colnames(ex.TCGA.t), function(x){
  unlist(paste(unlist(str_split(x,"\\."))[1:3],collapse = "-"))
}))
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\prognosis")
meta.2<-data.table::fread(file = "HCCDB15_time.txt",header = T)
rownames(meta.2)<-meta.2$Id
meta.2=data.frame(row.names = meta.2$Id,futime=meta.2$futime,fustat=meta.2$fustat)

d<-meta.2[a[which(duplicated(a))],]
rownames(d)<-paste0(a[which(duplicated(a))],c(1:3))
meta.2<-rbind(meta.2,d)
a[which(duplicated(a))]<-paste0(a[which(duplicated(a))],c(1:3))
setdiff(rownames(meta.2),a)
meta.2=meta.2[a,]
meta.2$fc<-TCGA.fc.4[colnames(ex.TCGA.t),]
colnames(ex.TCGA.t)=a
rownames(macro.TCGA)=a
rownames(cd8t.TCGA)=a


meta.t=meta[meta$tissue=="Tumor",]%>%as.data.frame()
rownames(meta.t)<-meta.t$icgc_specimen_id
meta.t<-meta.t[colnames(ex.ICGC.t),]

meta.3<-data.table::fread(file = "HCCDB18_time.txt",header = T)
rownames(meta.3)<-meta.3$id
meta.3=data.frame(row.names = meta.3$id,futime=meta.3$futime,fustat=meta.3$fustat)
index<-intersect(meta.t$icgc_specimen_id,rownames(ICGC.fc.4))
d<-meta.3[meta.t$icgc_donor_id,]
d1<-ICGC.fc.4[colnames(ex.ICGC.t),]
meta.t<-cbind(meta.t,d,d1)

ex.ICGC.t<-ex.ICGC.t[,meta.t$icgc_specimen_id]

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")

index<-intersect(Genes$Feature,rownames(ex.ICGC.t))
ex.ICGC.t<-ex.ICGC.t[index,]
index<-intersect(Genes$Feature,rownames(ex.TCGA.t))
ex.TCGA.t<-ex.TCGA.t[index,]

e<-ex.ICGC.t
name<-"ICGC.tumor"
set="ICGC"
group_list<-meta.2$fc
pcdim<-40

nc<-PCAnc(e,name,set,group_list,pcdim)

nc<-2
seed=1

df<-Dividgroup(e,name,set,group_list,pcdim,nc,seed = 1)

dfump<-df[,c(1:4,9:ncol(df))]
dftsne<-df[,5:ncol(df)]

cl<-as.factor(dfump[,5])
clname<-colnames(dfump)[5]
jpeg(file=paste0(name,"\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
UMAP(dfump,set,name,cl,clname,col=col4)
dev.off()
jpeg(file=paste0(name,"\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
TSNE(dftsne,set,name,cl,clname,col = col4)
dev.off()

for(i in 6:8){
  cl<-as.factor(dfump[,i])
  clname<-colnames(dfump)[i]
  jpeg(file=paste0(name,".\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  UMAP(dfump,set,name,cl,clname,col=cols)
  dev.off()
  jpeg(file=paste0(name,".\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  TSNE(dftsne,set,name,cl,clname,col = cols)
  dev.off()
}

for(i in 9:12){
  cl<-as.factor(dfump[,i])
  clname<-colnames(dfump)[i]
  jpeg(file=paste0(name,".\\log.",clname,".umap.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  UMAP(dfump,set,name,cl,clname,col=colsurv)
  dev.off()
  jpeg(file=paste0(name,".\\log.",clname,".tsne.",set,".jpeg"),width =6,height =5,units = "in", res = 1000)
  TSNE(dftsne,set,name,cl,clname,col = colsurv)
  dev.off()
}

#df<-cbind(dfump,macro.TCGA,cd8t.TCGA)
df<-cbind(df,macro.ICGC,cd8t.ICGC)
df$rank1=ifelse(df$CIBERSORT_M2>=median(df$CIBERSORT_M2),"High","Low")
df$rank2=ifelse(df$xCell_M2>=median(df$xCell_M2),"High","Low")
df$rank3=ifelse(df$CIBERSORT_M1>=median(df$CIBERSORT_M1),"High","Low")
df$rank4=ifelse(df$xCell_M1>=median(df$xCell_M1),"High","Low")
df$rank5=ifelse(df$MFP_M >=median(df$MFP_M),"High","Low")
df$rank6=ifelse(df$MFP_M1 >=median(df$MFP_M1),"High","Low")
table(df$rank1,df$c3)

cl=c()
for(i in 1:ncol(df)){
  a=class(df[,i])
  cl=c(cl,a)
}

a=chisq.test(df$fc,df$c3)
a$p.value
a=c(a$statistic,a$p.value)
chisq=data.frame(a)
for(i in colnames(df)[cl!="numeric"]){
  a=chisq.test(df$fc,df[,i])
  a=c(a$statistic,a$p.value)
  chisq=cbind(chisq,data.frame(a))
}
chisq=chisq[,-1]
colnames(chisq)=colnames(df)[cl!="numeric"]

library(ggalluvial)
data=df[,c("fc","rank3","rank4","rank1","rank2","rank5","rank6","c3","c2")]
data$fc=factor(data$fc,levels = c("D","F","IE/F","IE"))
data$rank4=factor(data$rank4,levels = c("Low","High"))
data$rank3=factor(data$rank3,levels = c("Low","High"))
data$rank1=factor(data$rank1,levels = c("Low","High"))
data$rank2=factor(data$rank2,levels = c("Low","High"))
levels(data$c3)<-c("CT","PT","AN")
levels(df$cluster3)<-c("CT","PT","AN")
unique(data$rank3)

if(sankey){
  jpeg(file=paste0(name,"\\feature.",set,".sankey.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data = data,
         aes(axis1=c3,axis2= fcad,label = after_stat(stratum))) +
    scale_x_discrete(limits = c("sample group", "MFP"), expand = c(.01, .05)) +
    geom_alluvium(aes(fill = c3)) +
    geom_stratum() + geom_text(stat = "stratum") +
    theme_minimal() +theme(panel.grid = element_blank())+
    scale_fill_manual(values =colsurv)
  dev.off()
  
  jpeg(file=paste0(name,"\\feature_",set,"_Proportiony.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data, aes(factor(c3)))+ 
    geom_bar(aes(fill = fcad), position = "fill")+ xlab("")+
    ylab("Proportion")+theme(legend.title=element_blank(),
                             panel.grid.major =element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.background = element_blank(),
                             panel.border = element_blank(),
                             strip.background.x = element_blank())+ 
    scale_fill_manual(values = c("#F98400","#E2D200","#00A08A","#5BBCD6"))+
    theme(axis.text.x=element_text(hjust=1))
  dev.off()
}

ex<-ex.ICGC.t#exprMatrix.3[rownames(ex.ICGC.t),]#ex.TCGA.t#all[rownames(ex.all.t),]
ex["CD74",]->cd74
ex["SPP1",]->spp1
ex["MIF",]->mif

arrydata<-rbind(MIF=mif,SPP1=spp1,CD74=cd74)%>%t()%>%as.data.frame()
arrydata<-log2(arrydata+1)%>%na.omit()%>%as.data.frame()
cor(arrydata[,1:2])

OPT<-function(x,data){
  ndata<-c(0,0,0)%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    s<-mean(data[i:(i+x),1])
    m<-mean(data[i:(i+x),2])
    cd<-mean(data[i:(i+x),3])
    col<-c(s,m,cd)%>%as.data.frame()
    ndata<-cbind(ndata,col)
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  #alist<-list()
  cdt<-c(0,0,0)%>%as.data.frame()
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-data.frame(c(cmat[1,2:3],cmat[2,3]))
    colnames(cc)<-i
    cdt<-cbind(cdt,cc)
    #alist=c(alist,list(ndt.i))
  }
  cdt<-cdt[,-1]%>%t()%>%as.data.frame()
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(4,arrydata)

arrydata<-cbind(t(ex[c("SPP1","MIF"),]),macro.ICGC,cd8t.ICGC)%>%as.data.frame()

OPT<-function(x,data){
  ndata<-c(rep(0,ncol(data)))%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    col<-c()
    for(j in 1:ncol(data)){
      s<-mean(data[i:(i+x),j])
      col<-c(col,s)
    }
    ndata<-cbind(ndata,as.data.frame(col))
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  cdt<-c(rep(0,ncol(data)*2-4))%>%as.data.frame()%>%t()
  colnames(cdt)<-c(paste("SPP1",colnames(data)[-2:-1],sep = "-"),paste("MIF",colnames(data)[-2:-1],sep = "-"))
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-t(data.frame(c(cmat[1,3:ncol(data)],cmat[2,3:ncol(data)])))
    colnames(cc)<-colnames(cdt)
    cdt<-rbind(cdt,cc)
  }
  cdt<-cdt[-1,]
  rownames(cdt)<-c(2:(y+1))
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(9,arrydata)

ce<-cor(t(ex))

l<-glm(arrydata.op[,1:2])
summary(l)

mylr = function(x,y){
  
  plot(x,y)
  
  x_mean = mean(x)
  y_mean = mean(y)
  xy_mean = mean(x*y)
  xx_mean = mean(x*x)
  yy_mean = mean(y*y)
  
  m = (x_mean*y_mean - xy_mean)/(x_mean^2 - xx_mean)
  b = y_mean - m*x_mean
  
  
  f = m*x+b# 线性回归方程
  
  lines(x,f)
  
  sst = sum((y-y_mean)^2)
  sse = sum((y-f)^2)
  ssr = sum((f-y_mean)^2)
  
  result = c(m,b,sst,sse,ssr)
  names(result) = c('m','b','sst','sse','ssr')
  
  return(result)
}

f<-mylr(arrydata.op[,2],arrydata.op[,1])
f['m']
f['b']
f['sse']+f['ssr']
f['sst']

R2 = f['ssr']/f['sst']

predict(l)
if(correlation){
  library(gglmannotate)
  jpeg(file=paste0(name,"\\seq-",set,"-MIF-SPP1-tumor.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data=arrydata.op, aes(x=MIF, y=SPP1)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (logTPM)") + ylab("SPP1 (logTPM)") + 
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
    geom_smooth(aes(x=MIF, y=SPP1),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                              place = "topleft")+
    annotate("text",x=5,y=8,label="R:0.64",parse=T,size=4)
  dev.off()
  
  bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
  jpeg(file=paste0(name,"\\cor.",set,".infi.jpeg"),width=8,height=6,units = "in", res = 1000)
  pheatmap::pheatmap(na.omit(arrycdt),#scale = "row",
                     #kmeans_k = 3,
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = F, cluster_col = T,
                     border_color ='white',
                     angle_col=45,
                     number_color="gray",
                     display_numbers =round(arrycdt, 2),
                     fontsize=5,
                     #fontsize_col=10,fontsize_row=5,
                     color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                               colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                     legend_breaks=seq(-1,1,0.2),
                     breaks=bk,
                     main ="correlation optimization") 
  dev.off()
}

if(circle){
  cora<-rbind(SPP1=arrycdt[9,1:(ncol(arrycdt)/2)],MIF=arrycdt[9,(ncol(arrycdt)/2+1):ncol(arrycdt)])%>%as.data.frame()
  colnames(cora)<-colnames(arrydata)[c(-2:-1)]
  set.seed(123)
  h<-pheatmap::pheatmap(t(cora),scale = "none",
                        show_colnames =T,
                        show_rownames = T,
                        cluster_row = T, cluster_col = F,
                        border_color =NA,
                        number_color="black",
                        display_numbers =round(t(cora), 2),
                        cellwidth = 30, cellheight = 30,
                        angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                        na_col = "grey",
                        color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                                  colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                        legend_breaks=seq(-1,1,0.5),
                        breaks=bk,
                        main ="correlation") 
  jpeg(file=paste0(name,"\\cor.",set,".jpeg"),width=6,height=8,units = "in", res = 1000)
  h
  dev.off()
  
  fdf<-t(cora)
  a<-c(18,6,13,12,14,17,9,15,8,11,4,3,5,10,16,7,1,2)
  name<-rownames(fdf)[h$tree_row$order]#
  fdf<-fdf[a,]
  
  library(reshape2)
  cirmat<-melt(fdf)
  head(cirmat)
  
  resn <- cirmat %>% filter(Var2 == 'MIF')
  rese <- cirmat %>% filter(Var2 == 'SPP1')
  
  resn$ang <- seq(from = (360/nrow(resn)) / 1.5,
                  to = (1.5* (360/nrow(resn))) - 360,
                  length.out = nrow(resn)) + 80
  
  resn$hjust <- 0
  resn$hjust[which(resn$ang < -90)] <- 1
  resn$ang[which(resn$ang < -90)] <- (180+resn$ang)[which(resn$ang < -90)]
  rese$hjust<-resn$hjust
  rese$ang<-resn$ang
  cirmat$var <- rep(1:nrow(resn),2)
  
  s_anno <- data.frame(value = c("SPP1","MIF"))
  s_anno$value <-factor(s_anno$value,levels = c("SPP1","MIF"))
  cirmat$Var2<-factor(cirmat$Var2,levels = c("SPP1","MIF"))
  
  library(gplots)
  library(ggdendro)
  library(ggh4x)
  library(ggnewscale)
  library(ggrepel)
  p1 <- ggplot() +
    geom_bar(data = s_anno,stat = 'identity',
             aes(x = 0,y = 1,fill = value),
             width = 1,
             color = NA) +
    scale_fill_manual(name = 'Gene',
                      values = c("MIF"="#E2D200", 
                                 "SPP1"="#F98400"))
  p2<-p1+ new_scale("fill")+
    geom_tile(data = cirmat[which(cirmat$Var2 == 'MIF'),],
              aes(x = 1:nrow(resn),y = 0.5,fill = value),
              color = 'white') +
    geom_tile(data = cirmat[which(cirmat$Var2 == 'SPP1'),],
              aes(x = 1:nrow(resn),y = 1.5,fill = value),
              color = 'white')+  scale_fill_gradient2(midpoint = 0,
                                                      low = "#5BBCD6",
                                                      mid = "white",
                                                      high = "#B40F20") +
    ylim(-1,5)+xlim(-1,20)
  
  name<-"ICGC.tumor"
  jpeg(file=paste0(name,"\\cor.",set,".circle-1.jpeg"),width=9,height=6,units = "in", res = 1000)
  jpeg(file="cor.circle-1.jpeg",width=9,height=6,units = "in", res = 1000)
  p2 + coord_polar(theta = 'x') +theme_void() +new_scale("colour")+
    geom_text(data = resn,
              aes(x = as.numeric(rownames(resn)),
                  y = 2.2,
                  label = as.character(Var1), angle = ang, hjust = hjust),
              size = 4)
    # geom_text(data = resn,
    #           aes(x = as.numeric(rownames(resn)),
    #               y = c(0.25),
    #               label = round(value,2), angle = ang, hjust = hjust),
    #           size = 3.5)+
    # geom_text(data = rese,
    #           aes(x = as.numeric(rownames(rese)),
    #               y = c(1.25),
    #               label = round(value,2),angle = ang, hjust = hjust),
    #           size = 4.5)#+NoLegend()
  # annotate("text", x = as.numeric(rownames(rese)), 
  #          y = 0.25, label = round(rese$value,2), color = "white", size = 4)
  dev.off()
}

df.infi<-cbind(macro.TCGA,cd8t.TCGA)
###################array-all#######
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
load("intergated_array.Rdata")
ann_col<-data.frame(sample=colnames(data),tissue=tissue,batch=batch)
rownames(ann_col)<-ann_col$sample
ann_col_t<-ann_col[ann_col$tissue=="Tumor",]

rm(data)
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge-N")
all_xcell=read.csv(".\\results\\all_xcell.csv", row.names = 1)
all_CIBERSORT=data.table::fread(".\\results\\all-CIBERSORT-Results.txt")
all_CIBERSORT<-CIBERSORT_NAME(all_CIBERSORT)
all_MCP<-read.csv(".\\results\\all_MCP.csv", row.names = 1)
all_xcell.t<-all_xcell[,ann_col$sample[ann_col$tissue=="Tumor"]]
all_CIBERSORT.t<-all_CIBERSORT[,ann_col$sample[ann_col$tissue=="Tumor"]]
all_MCP.t<-all_MCP[,ann_col$sample[ann_col$tissue=="Tumor"]]
array.fc.4<-read.csv(".\\results\\array_tumor_log_final_clusters_4.csv",row.names = 1)

array.mfp<-read.csv(".\\results\\array_all_log_signature_scores_scaled.csv",row.names = 1)
array.mfp<-array.mfp[rownames(ann_col),]
array.mfp.t<-array.mfp[ann_col$sample[ann_col$tissue=="Tumor"],]
array.mfp.add<-read.csv(".\\results\\array_all_log_signature_scores_scaled_add.csv",row.names = 1)
array.mfp.add<-array.mfp.add[rownames(ann_col),]
array.mfp.add.t<-array.mfp.add[ann_col$sample[ann_col$tissue=="Tumor"],]
array.fc.4.add<-read.csv(".\\results\\array_tumor_log_final_clusters_4_add.csv",row.names = 1)
array.fc.4<-array.fc.4[ann_col_t$sample,]
array.fc.4.add<-array.fc.4.add[ann_col_t$sample,]

ann_col_t$fc4ad<-array.fc.4.add
ann_col_t$fc4<-array.fc.4

index=intersect(rownames(all),Genes$Feature)
ex.all=all[index,]
ex_z.all=Zscore(all[index,])
ex.all.t<-ex.all[,ann_col_t$sample]
ex_z.t<-ex_z.all[,ann_col_t$sample]

e<-ex.all.t
name<-"all.tumor"
set="array"
group_list<-ann_col_t$fc4ad
pcdim<-40

#nc<-PCAnc(e,name,set,group_list,pcdim)

nc<-2
seed=1
df<-Dividgroup(e,name,set,group_list,pcdim,nc,seed)

cluster3.w<-CLUSTER3.wt(ex.all.t)

cl=c()
for(i in 1:ncol(df)){
  a=class(df[,i])
  cl=c(cl,a)
}

a=chisq.test(df$fc,df$c3)
a$p.value
a=c(a$statistic,a$p.value)
chisq=data.frame(a)
for(i in colnames(df)[cl!="numeric"]){
  a=chisq.test(df$c3,df[,i])
  a=c(a$statistic,a$p.value)
  chisq=cbind(chisq,data.frame(a))
}
chisq=chisq[,-1]
colnames(chisq)=colnames(df)[cl!="numeric"]

library(ggalluvial)
data=df[,c("fc","fcad","rank3","rank4","rank1","rank2","rank5","rank6","c3","c2")]
data$fc=factor(data$fc,levels = c("D","F","IE/F","IE"))
data$fcad=factor(data$fcad,levels = c("D","F","IE/F","IE"))
data$rank4=factor(data$rank4,levels = c("Low","High"))
data$rank3=factor(data$rank3,levels = c("Low","High"))
data$rank1=factor(data$rank1,levels = c("Low","High"))
data$rank2=factor(data$rank2,levels = c("Low","High"))
levels(data$c3)<-c("CT","PT","AN")
levels(df$cluster3)<-c("CT","PT","AN")
unique(data$rank3)

ex<-ex.all.t#all[rownames(ex.all.t),]
ex["CD74",]->cd74
ex["SPP1",]->spp1
ex["MIF",]->mif

arrydata<-rbind(MIF=mif,SPP1=spp1,CD74=cd74)%>%t()%>%as.data.frame()
#arrydata<-log2(arrydata+1)%>%na.omit()%>%as.data.frame()
cor(arrydata[,1:2])

OPT<-function(x,data){
  ndata<-c(0,0,0)%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    s<-mean(data[i:(i+x),1])
    m<-mean(data[i:(i+x),2])
    cd<-mean(data[i:(i+x),3])
    col<-c(s,m,cd)%>%as.data.frame()
    ndata<-cbind(ndata,col)
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  #alist<-list()
  cdt<-c(0,0,0)%>%as.data.frame()
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-data.frame(c(cmat[1,2:3],cmat[2,3]))
    colnames(cc)<-i
    cdt<-cbind(cdt,cc)
    #alist=c(alist,list(ndt.i))
  }
  cdt<-cdt[,-1]%>%t()%>%as.data.frame()
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(4,arrydata)

arrydata<-cbind(data.frame(SPP1=spp1,MIF=mif),macro,cd8t)%>%as.data.frame()

OPT<-function(x,data){
  ndata<-c(rep(0,ncol(data)))%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    col<-c()
    for(j in 1:ncol(data)){
      s<-mean(data[i:(i+x),j])
      col<-c(col,s)
    }
    ndata<-cbind(ndata,as.data.frame(col))
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  cdt<-c(rep(0,ncol(data)*2-4))%>%as.data.frame()%>%t()
  colnames(cdt)<-c(paste("SPP1",colnames(data)[-2:-1],sep = "-"),paste("MIF",colnames(data)[-2:-1],sep = "-"))
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-t(data.frame(c(cmat[1,3:ncol(data)],cmat[2,3:ncol(data)])))
    colnames(cc)<-colnames(cdt)
    cdt<-rbind(cdt,cc)
  }
  cdt<-cdt[-1,]
  rownames(cdt)<-c(2:(y+1))
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(9,arrydata)

ce<-cor(t(ex))

if(correlation){
  library(gglmannotate)
  jpeg(file=paste0(name,"\\array-MIF-SPP1-tumor.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data=scdata.op, aes(x=MIF, y=SPP1)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (TPM)") + ylab("SPP1 (TPM)") + 
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
    geom_smooth(aes(x=MIF, y=SPP1),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                             place = "topleft")+
    annotate("text",x=0,y=75,label="R:0.50",parse=T,size=4)
  dev.off()
  
  jpeg(file=paste0(name,"\\cor.all.intergrated.infi.jpeg"),width=8,height=6,units = "in", res = 1000)
  pheatmap::pheatmap(na.omit(arrycdt),#scale = "row",
                     #kmeans_k = 3,
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = F, cluster_col = T,
                     border_color ='white',
                     angle_col=45,
                     number_color="gray",
                     display_numbers =round(arrycdt, 2),
                     fontsize=5,
                     #fontsize_col=10,fontsize_row=5,
                     color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                               colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                     legend_breaks=seq(-1,1,0.2),
                     breaks=bk,
                     main ="correlation optimization") 
  dev.off()
}
save(list = ls(),file = "array-intergrated.Rdata")
###################array-divide#######
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge-N")
ex.1.t<-ex.all.t[,rownames(time.1)]
ex.2.t<-ex.all.t[,rownames(time.2)]
name<-"1.tumor"
set<-"1"
cluster3.w<-CLUSTER3.wt(ex.1.t)
cluster3.w<-CLUSTER3.wt(ex.2.t)
##################sankey#################
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")
ICGC.cli<-data.table::fread("ICGC-clinical.txt")%>%as.data.frame()
rownames(ICGC.cli)<-ICGC.cli$icgc_specimen_id
ICGC.cli<-ICGC.cli[meta.t$icgc_specimen_id,]
table(ICGC.cli$tumour_grade)
table(ICGC.cli$tumour_stage)
ICGC.cli$tumour_grade<-ifelse(ICGC.cli$tumour_grade%in%c("I","I-II","I-III"),'G1',
                              ifelse(ICGC.cli$tumour_grade%in%c("II","II-I"),"G2",
                                     ifelse(ICGC.cli$tumour_grade%in%c("II-III","III"),"G3",
                                            ifelse(ICGC.cli$tumour_grade=="IV","G4","unknow"))))

TCGA.cli<-data.table::fread("TCGA-clinical.txt")%>%as.data.frame()
rownames(TCGA.cli)<-TCGA.cli$Id
TCGA.cli<-TCGA.cli[rownames(meta.2),]
table(TCGA.cli$stage)
table(TCGA.cli$grade)
TCGA.cli$stage<-ifelse(TCGA.cli$stage%in%c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC"),3,
                       ifelse(TCGA.cli$stage%in%c("Stage IV","Stage IVA","Stage IVB"),4,
                              ifelse(TCGA.cli$stage=="Stage II",2,
                                     ifelse(TCGA.cli$stage=="Stage I",1,"unknow"))))

cli<-data.frame(row.names = colnames(ex.all.t),grade=c(TCGA.cli$grade,ICGC.cli$tumour_grade),
                stage=c(TCGA.cli$stage,ICGC.cli$tumour_stage))
cli$cluster<-cluster3.w
table(cluster3.w,cli$stage)

GSE14520<-data.table::fread("GSE14520_Extra_Supplement.txt")%>%as.data.frame()
rownames(GSE14520)<-GSE14520$Affy_GSM
GSE14520<-GSE14520[rownames(GSE14520)%in%colnames(ex.all.t),]
table(GSE14520$`TNM staging`)
GSE14520$`TNM staging`<-ifelse(GSE14520$`TNM staging`%in%c("IIIA","IIIB", "IIIC"),3,
                               ifelse(GSE14520$`TNM staging`=="II",2,
                                      ifelse(GSE14520$`TNM staging`=="I",1,"unknow")))
cli.a<-data.frame(cluster=cluster3.w[rownames(GSE14520)],stage=GSE14520$`TNM staging`,row.names = rownames(GSE14520))


cli$grade[which(is.na(cli$grade))]="unknow"
cli$stage[which(is.na(cli$stage))]="unknow"
cli$grade<-plyr::mapvalues(x=cli$grade,from =c("G1","G2","G3","G4","unknow"),to=c("1","2","3","4","unknow") )
cli$grade<-factor(cli$grade,levels=c("unknow","1","2","3","4"))
cli$stage<-factor(cli$stage,levels=c("unknow","1","2","3","4"))
cli$grade2<-ifelse(cli$grade%in%c(1,2),"G1-G2",
                   ifelse(cli$grade%in%c(3,4),"G3-G4","unknown"))
cli$stage2<-ifelse(cli$stage%in%c(1,2),"I-II",
                   ifelse(cli$stage%in%c(3,4),"III-IV","unknown"))

cli.a$stage<-factor(cli.a$stage,levels=c("unknow","1","2","3","4"))
cli.a$stage<-plyr::mapvalues(x=cli.a$stage,from=c("unknow","1","2","3","4"),to=c("unknow","I","II","III","IV"))

jpeg(file="feature_TCGA_stage_grade_sankey-1.jpeg",width=6,height=4,units = "in", res = 600)
ggplot(data = cli,
       aes(axis1=cluster,axis2=grade2,axis3= stage2,label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sample group", "Edmondson grade","TNM stage"), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = cluster)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +theme(panel.grid = element_blank())+
  scale_fill_manual(values =colsurv)
dev.off()

jpeg(file="feature_TCGA_grade_sankey.jpeg",width=4,height=4,units = "in", res = 600)
ggplot(data = cli.T,
       aes(axis1=cluster,axis2=grade,label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sample group", "Edmondson grade"), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = cluster)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +theme(panel.grid = element_blank())+
  scale_fill_manual(values =colsurv)
dev.off()

jpeg(file="feature_TCGA_stage_sankey.jpeg",width=4,height=4,units = "in", res = 600)
ggplot(data = cli.T,
       aes(axis1=cluster,axis2=stage,label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sample group", "TNM stage"), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = cluster)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +theme(panel.grid = element_blank())+
  scale_fill_manual(values =colsurv)
dev.off()

ncli<-data.frame(row.names=c(rownames(cli.a),rownames(cli)),
                 stage=c(cli.a$stage,cli$stage),cluster=c(as.character(cli.a$cluster),as.character(cli$cluster)))

jpeg(file="TCGA_grade_Proportion.jpeg",width=4,height=4,units = "in", res = 2000)
ggplot(cli.T, aes(factor(cluster)))+ 
  geom_bar(aes(fill = stage), position = "fill")+ xlab("")+
  ylab("Proportion")+theme(legend.title=element_blank(),
                           panel.grid.major =element_blank(), 
                           panel.grid.minor = element_blank(), 
                           panel.background = element_blank(),
                           panel.border = element_blank(),
                           strip.background.x = element_blank())+ 
  scale_fill_manual(values = col5)+
  theme(axis.text.x=element_text(hjust=1))
dev.off()

jpeg(file="GSE14520_stage_Proportion-1.jpeg",width=4,height=4,units = "in", res = 2000)
ggplot(cli.a, aes(factor(stage)))+ 
  geom_bar(aes(fill = cluster), position = "fill")+ xlab("")+
  ylab("Proportion")+theme(legend.title=element_blank(),
                           panel.grid.major =element_blank(), 
                           panel.grid.minor = element_blank(), 
                           panel.background = element_blank(),
                           panel.border = element_blank(),
                           strip.background.x = element_blank())+ 
  scale_fill_manual(values = colsurv)+
  theme(axis.text.x=element_text(hjust=1))
dev.off()
cli.T<-cli[rownames(meta.2),]
cli.I<-cli[rownames(meta.t),]

annotation_col = data.frame(
  #km=factor(df$km),
  #MFP=factor(df$fc),
  #batch=factor(ann_col_t$batch)
)
c("I"="#D3D3D3","II"="#C89197","III"="#BE505B","unknow"="white")
anncols<-list(stage=col5.1,cluster=colsurv)
rownames(annotation_col) =colnames(ex.all.t)

table(dfump$km,dftsne$km)

h<-pheatmap::pheatmap(t(ex.all.t),scale = "column",
                      #kmeans_k = 3,
                      show_colnames =F,
                      show_rownames = F,
                      #annotation_col = cli.a,
                      annotation_row  =cli,
                      cluster_row = T, 
                      cluster_col = T,
                      annotation_colors = anncols,
                      # border_color ='white',
                      # angle_col=45,#fontsize=10,fontsize_col=10,fontsize_row=6,
                      color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdYlBu")))(100),
                      # legend = F,
                      main ="score")
jpeg(file="seq_heatmap.jpeg",width=4,height=4,units = "in", res = 2000)
h
dev.off()
ex.all.t<-ex.all.t[,rownames(cli.a)]
################survive####################
df$batch<-batch.t

library(survival)
library(survivalROC)
library(beeswarm)
library("survminer")
library(preprocessCore)
meta.2$cluster<-df$c3[df$batch=="TCGA"]%>%as.factor()
meta.t$cluster<-df$c3[df$batch=="ICGA"]%>%as.factor()

meta.2$cluster<-meta.2$fc%>%as.factor()
meta.t$cluster<-meta.t$d1%>%as.factor()
meta.t$cluster<-df$c3%>%as.factor()
meta.2$cluster<-df$c3


rt<-meta.2
rt<-na.omit(rt)
rt$futime=rt$futime/30
table(rt$cluster)
max(rt$futime)

#rt<-rt[rt$futime<=36,]
rt.t<-rt[rt$cluster!="CT",]

fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
diff=survdiff(Surv(futime, fustat) ~cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=length(diff$n)-1)
pValue
diff
# pValue=signif(pValue,4)

#0.005053787
#0.0004279965
levels(df$c3)
colsurv=c("CT"="#FF0000","AN"="#00A08A","PT"="#F98400")
jpeg(file=paste0(name,"\\0907-survival-",set,".jpeg"),width =5,height =6,units = "in", res = 2000)
ggsurvplot(fit, 
           data=rt,
           #conf.int=TRUE,
           palette=c("#00A08A","#F98400","#FF0000"),#c("#F98400","#E2D200","#00A08A","#5BBCD6"),#
           pval=paste0("p=",pValue),
           pval.size=4,
           #legend.labs=c("High", "Low"),
           #legend.title="Expression",
           break.time.by = 10,
           risk.table.title="",
           risk.table=TRUE,
           risk.table.height=.30,
           #xlim=c(0,24),
           ylim=c(0,1),
           xlab="Time(months)")
dev.off()

save(list = ls(),file = "D:\\Bioinfrolf\\Bioinfrolf\\tmpRdata/TCGA_ICGC_seperated.Rdata")

ex<-t(ex.all.t[,batch.t=="ICGA"])
ex<-t(ex.TCGA.t)
meta<-cbind(meta.2,ex[,c("MIF","SPP1")])

df1<-meta
max(df1$futime)

df1$futime=df1$futime/30
df1<-df1[order(df1$MIF,decreasing = T),]
df1$levels<-ifelse(df1$MIF>=median(df1$MIF),"High expression","Low expression")

#df1<-df1[df1$futime<=60,]
outTab=data.frame()
for (i in 10:200){
  for(j in 10:200){
    rt<-df1[c(1:i,(nrow(df1)-j+1):nrow(df1)),]
    fit <- survfit(Surv(futime, fustat) ~ levels, data = rt)
    diff=survdiff(Surv(futime, fustat) ~levels,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    pValue=signif(pValue,4)
    out=data.frame(h=i,l=j,p=pValue)
    outTab=rbind(outTab,out)
  }
}

outTab<-outTab[order(outTab$p,decreasing = F),]
outTab$d<-abs(outTab$h-outTab$l)
outTab$a<-outTab$h+outTab$l
outTab<-outTab[outTab$p<0.05,]
outTab.1<-outTab[outTab$d<50,]
outTab.1<-outTab.1[ outTab.1$a >=150,]
outTab.1<-outTab.1[order(outTab.1$p,decreasing = F),]
write.csv(outTab,"MIF-pvalue_TCGA.csv")

i=73
j=112
rt<-df1[c(1:i,(nrow(df1)-j+1):nrow(df1)),]
fit <- survfit(Surv(futime, fustat) ~ levels, data = rt)
diff=survdiff(Surv(futime, fustat) ~levels,data = rt)
pValue=1-pchisq(diff$chisq,df=length(diff$n)-1)
pValue=signif(pValue,4)

jpeg(file=paste0(name,"\\survival-TCGA-MIF.jpeg"),width =5,height =6,units = "in", res = 1000)
ggsurvplot(fit, 
           data=rt,
           #conf.int=TRUE,
           palette=c("#B40F20","#5BBCD6"),
           pval=paste0("p=",pValue),
           pval.size=4,
           legend.labs=c("High", "Low"),
           legend.title="Expression",
           break.time.by = 6,
           risk.table.title="",
           risk.table=TRUE,
           risk.table.height=.30,
           #xlim=c(0,24),
           #ylim=c(0.2,1),
           xlab="Time(months)")
dev.off()

save(list = ls(),file = "TCGA.Rdata")
#################DEG###################################

all<-ex.all#na.omit(exprMatrix.3)
conNum=50           
treatNum=374
index<-intersect(Genes$Feature,rownames(all))
barplot(1:8,col=col8)
data=all[index,]
dim(data)

patient<-colnames(data)
table(tissue)
#which(tissue=="Normal")
ajn<-patient[tissue=="Normal"]
te<-patient[tissue=="Tumor"]

data=cbind(data[,ajn],data[,te])%>%as.matrix()
treatNum=length(te)
conNum=length(ajn)  
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))

i="SPP1"
for(i in row.names(data)){
  #geneName=unlist(strsplit(i,"\\|",))[1]
  #geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  t(rt)
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(as.numeric(data[i,1:conNum]))
  treatGeneMeans=mean(as.numeric(data[i,(conNum+1):ncol(data)]))
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(as.numeric(data[i,1:conNum]))
  treatMed=median(as.numeric(data[i,(conNum+1):ncol(data)]))
  diffMed=treatMed-conMed
  logFC.med=log2(treatMed)-log2(conMed)
  #if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){   } 
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,
                              conMed=conMed,treatMed=treatMed,logFC.med=logFC.med,pValue=pvalue))

}

pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
outTab$pValue<-as.numeric(outTab$pValue)
outTab=outTab[order(outTab$pValue,decreasing = F),]
outTab[which(outTab$gene%in%c("MIF","CD74","SPP1")),]
# outTab=outTab[outTab$pValue<0.05,]
# outTab=outTab[outTab$logFC>0,]
write.table(outTab,file="T&N-DE-array-all.txt",row.names=F,quote=F)
outTab.T<-outTab
save(outTab,grade,data,file="T&N-DE-array.Rdata")

library(limma)
library(umap)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
which(c("EGF","NOX4")%in%G2)
rm(list = ls())
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\seq")
sDEs=read.table("T&N-DE-seq-all.txt",header = T,row.names = 1)
sDEs[c("MIF","CD74","SPP1"),]

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge-N")
sDEa=read.table("T&N-DE-array-all.txt",header = T,row.names = 1)
sDEa[c("MIF","CD74","SPP1"),]

load("D:/Bioinfrolf/Bioinfrolf/valid/HCCDB/scGene.Rdata")
sDEs$p=ifelse(sDEs$fdr<0.05,1,0)
sDEa$p=ifelse(sDEa$fdr<0.05,1,0)

int=intersect(rownames(sDEs),rownames(sDEa))
dimnames=list(int,c("seq","array"))
intdf=matrix(rep(0,length(int)*2),nrow=length(int),dimnames=dimnames)%>%as.data.frame()
intdf$seq=sDEs[int,"logFC"]
intdf$array=sDEa[int,"logFC"]

intersect(rownames(sDEs)[sDEs$p==1],rownames(sDEa)[sDEa$p==1])
intdf[which(rownames(intdf)%in%c("MIF","CD74","SPP1")),]

intdf$result = as.factor(ifelse(intdf$seq>0.5 & intdf$array > 0.2,"UP",
                                ifelse(intdf$seq<(-0.5) & intdf$array<(-0.2),
                                       "DOWN",'NOT')))
intdf$Label = c(rep(NA,nrow(intdf)))
intdf$Label[which(rownames(intdf)%in%c("MIF","SPP1"))]=c("SPP1","MIF")
co=cor(as.numeric(intdf$seq),as.numeric(intdf$array))


hr<-list(t=rownames(intdf)[intdf$result=="UP"],
         g1=intersect(rownames(intdf),G1),
         g2=intersect(rownames(intdf),G2))
intersect(hr$t,hr$g2)%>%length()

lr<-list(t=rownames(intdf)[intdf$result=="DOWN"],
         g=intersect(rownames(intdf),G0))
intersect(lr$t,lr$g)%>%length()
length(lr$t)
sD
library(VennDiagram)   
library(gplots)
col3<-c("#B40F20","#FF0000","#F98400")
col2=c("#00A08A","#5BBCD6")
col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6",
         "#E2D200","#B40F20","#273046","#FD6467")
barplot(1:8,col=col8)
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB")
library(eulerr)
vd <- euler(c(A= (length(hr$t)-19-8),
              B=(length(hr$g2)-19),
              C=(length(hr$g1)-8), 
              "A&B" = 19, "A&C" =8, "B&C" = 0,
              "A&B&C" = 0))
length(hr$g1)+length(hr$g2)+length(lr$g)
pvalue <- phyper(27-1,178,341-178,34,lower.tail=F)
qvalue <- p.adjust(pvalue,method='fdr')
qvalue
phyper(27-1,178,341-178,31,lower.tail = F)%>%p.adjust()#2.212886e-05

jpeg(file="VN-h-2.25.jpeg",width=3,height=3,units = "in", res = 1000)
plot(vd,
     fills = list(fill = col3, alpha = 0.4),
     labels = list(col = "white", font = 1,size=0.1,family="serif"), 
     edges = T,family="serif",
     quantities = TRUE)
dev.off()

vd <- euler(c(A= (length(lr$g)-14),
              B=(length(lr$t)-14),
              "A&B" = 14))
phyper(14-1,length(lr$g),341-163,29,lower.tail = F)
phyper(18-1,length(lr$g),341-163,37,lower.tail = F)#0.525025

jpeg(file="VN-l-2.25.jpeg",width=3,height=3,units = "in", res = 1000)
plot(vd,
     fills = list(fill = col2, alpha = 0.4),
     labels = list(col = "white", font = 1), 
     edges = T,family = "serif",
     quantities = TRUE)
dev.off()

library(gglmannotate)
jpeg(file="D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\Tumor-DE-means-1.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(intdf,aes(x=array,y=seq))+geom_point(aes(color = result),size=2)+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text.x =element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),#,family = "serif"
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12)
  )+xlab('logFC(array)') + ylab('logFC(seq)')+ggtitle("DEG")+
  scale_color_manual(values =c("#5BBCD6","lightgray","#B40F20"))+
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed",size=0.5)+ 
  geom_vline(xintercept= c(-0.2,0.2),linetype="dashed",size=0.5)+
  geom_label_repel(aes(label=Label),size=5,color="black",# fontface="bold",
                   box.padding=unit(1, "lines"), point.padding=unit(0.3, "lines"),
                   segment.colour = "black",na.rm = TRUE,show.legend = FALSE,force = T,
                   nudge_x = 0.1, nudge_y = 0.6)+
  geom_smooth(aes(x=array, y=seq),method = lm,linetype="dashed",se=T,size=0.5, formula=y ~ x,colour="black",
              fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                           place = "topleft",min.size=4)
dev.off()

#################array#############################
rm(list = ls())
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB6_mRNA")
exprMatrix.1 = read.table(file = "HCCDB6_GSE14520-GPL3921_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
meta.1<-data.table::fread(file = "HCCDB6.sample.txt",header = F)
meta.1<-t(meta.1)%>%as.data.frame()
colnames(meta.1)<-meta.1[1,]
meta.1<-meta.1[-1,]
rownames(meta.1)<-meta.1$GEO_ID
meta.1<-meta.1[order(meta.1$TYPE,decreasing = T),]
max(which(meta.1$TYPE=="HCC"))
meta.1.t<-meta.1[meta.1$TYPE=="HCC",]
time.1<-read.table("HCCDB6.patient.txt",header = T,row.names = 1)
time.1<-time.1%>%t()%>%as.data.frame()
time.1<-time.1[gsub("-","\\.",meta.1.t$PATIENT_ID),]
rownames(time.1)<-rownames(meta.1.t)

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB17_mRNA")
exprMatrix.2 = read.table(file = "HCCDB17_GSE76427_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
meta.2<-data.table::fread(file = "HCCDB17.sample.txt",header = F)
meta.2<-t(meta.2)%>%as.data.frame()
colnames(meta.2)<-meta.2[1,]
meta.2<-meta.2[-1,]
rownames(meta.2)<-meta.2$GEO_ID
meta.2<-meta.2[order(meta.2$TYPE,decreasing = T),]
max(which(meta.2$TYPE=="HCC"))
meta.2.t<-meta.2[meta.2$TYPE=="HCC",]
time.2<-read.table("HCCDB17.patient.txt",header = T,row.names = 1)
time.2<-time.2%>%t()%>%as.data.frame()
rownames(time.2)<-time.2$PATIENT_ID
time.2<-time.2[meta.2.t$PATIENT_ID,]
rownames(time.2)<-rownames(meta.2.t)

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB7_mRNA\\GSE10141-HCC-array")
exprMatrix.3 = read.table(file = "GSE10141_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
time.3<-read.table("GSE10141_time.txt",header = T)
rownames(time.3)<-time.3$id
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\data\\HCCDB7_mRNA")
meta.3<-data.table::fread(file = "HCCDB7.sample.txt",header = F)
meta.3<-t(meta.3)%>%as.data.frame()
colnames(meta.3)<-meta.3[1,]
meta.3<-meta.3[-1,]
rownames(meta.3)<-meta.3$GEO_ID
meta.3<-meta.3[order(meta.3$TYPE,decreasing = T),]
max(which(meta.3$TYPE=="HCC"))
meta.3.t<-meta.2[rownames(time.3),]

ex.1.t<-exprMatrix.1[,rownames(meta.1.t)]
ex.2.t<-exprMatrix.2[,rownames(meta.2.t)]
ex.3.t<-exprMatrix.3[,rownames(time.3)]

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge-N")

#write.csv(t(log2(ex.3.t+1)),file = ".\\data\\GSE10141_tumor_symbol_logtpm.csv")
# write.csv(t(log2(ex.1.t+1)),file = ".\\data\\GSE14520_tumor_symbol_logtpm.csv")
# write.csv(t(log2(ex.2.t+1)),file = ".\\data\\GSE76427_tumor_symbol_logtpm.csv")
# write.csv(t(ex.1.t),file = ".\\data\\GSE14520_tumor_symbol_tpm.csv")
# write.csv(t(ex.2.t),file = ".\\data\\GSE76427_tumor_symbol_tpm.csv")
# rt=cbind(ID=rownames(ex.1.t),as.data.frame(ex.1.t))
# write.table(rt,file = "GSE14520_tumor_symbol_tpm.txt",quote = F,row.names = F)
# rt=cbind(ID=rownames(ex.2.t),as.data.frame(ex.2.t))
# write.table(rt,file = "GSE76427_tumor_symbol_tpm.txt",quote = F,row.names = F)


# xcell.e1=xCellAnalysis(ex.1.t)
# write.csv(xcell.e1,file = ".\\results\\GSE14520_xcell.csv",quote = F)
# e1_MCP<- MCPcounter.estimate(ex.1.t,featuresType=c("HUGO_symbols")[1],
#                                     probesets=probesets,
#                                     genes=genes
# )
# write.csv(e1_MCP,file = ".\\results\\GSE14520_MCP.csv",quote = F)
# results0=CIBERSORT(".\\data\\ref.txt", ".\\GSE14520_tumor_symbol_tpm.txt", perm=1000,QN=TRUE,name="GSE14520_")
e1_xcell=read.csv(".\\results\\GSE14520_xcell.csv", row.names = 1)
e1_MCP<-read.csv(".\\results\\GSE14520_MCP.csv", row.names = 1)
e1_CIBERSORT=data.table::fread(".\\results\\GSE14520_CIBERSORT-Results.txt")
e1_CIBERSORT<-CIBERSORT_NAME(e1_CIBERSORT)
e1.mfp<-read.csv(".\\results\\GSE14520_log_signature_scores_scaled.csv",row.names = 1)
e1.mfp<-e1.mfp[rownames(meta.1.t),]
# xcell.e2=xCellAnalysis(ex.2.t)
# write.csv(xcell.e2,file = ".\\results\\GSE76427_xcell.csv",quote = F)
# e2_MCP<- MCPcounter.estimate(ex.2.t,featuresType=c("HUGO_symbols")[1],
#                                     probesets=probesets,
#                                     genes=genes
# )
# write.csv(e2_MCP,file = ".\\results\\GSE76427_MCP.csv",quote = F)
# results0=CIBERSORT(".\\data\\ref.txt", ".\\GSE76427_tumor_symbol_tpm.txt", perm=1000,QN=TRUE,name="GSE76427_")
e2_xcell=read.csv(".\\results\\GSE76427_xcell.csv", row.names = 1)
e2_MCP<-read.csv(".\\results\\GSE76427_MCP.csv", row.names = 1)
e2_CIBERSORT=data.table::fread(".\\results\\GSE76427_CIBERSORT-Results.txt")
e2_CIBERSORT<-CIBERSORT_NAME(e2_CIBERSORT)
e2.mfp<-read.csv(".\\results\\GSE76427_log_signature_scores_scaled.csv",row.names = 1)
e2.mfp<-e2.mfp[rownames(meta.2.t),]

# xcell.e3=xCellAnalysis(ex.3.t)
# write.csv(xcell.e3,file = ".\\results\\GSE10141_xcell.csv",quote = F)
# e3_MCP<- MCPcounter.estimate(ex.3.t,featuresType=c("HUGO_symbols")[1],
#                                     probesets=probesets,
#                                     genes=genes
# )
# write.csv(e3_MCP,file = ".\\results\\GSE10141_MCP.csv",quote = F)
# results0=CIBERSORT(".\\data\\ref.txt", ".\\GSE10141_symbol_tpm.txt", perm=1000,QN=TRUE,name="GSE10141-")
e3_xcell=read.csv(".\\results\\GSE10141_xcell.csv", row.names = 1)
e3_MCP<-read.csv(".\\results\\GSE10141_MCP.csv", row.names = 1)
e3_CIBERSORT=data.table::fread(".\\results\\GSE10141-CIBERSORT-Results.txt")
e3_CIBERSORT<-CIBERSORT_NAME(e3_CIBERSORT)
e3.mfp<-read.csv(".\\results\\GSE10141_log_signature_scores_scaled.csv",row.names = 1)
e3.mfp<-e3.mfp[rownames(time.3),]

index<-intersect(Genes$Feature,rownames(ex.1.t))
ex.1.t<-ex.1.t[index,]
index<-intersect(Genes$Feature,rownames(ex.2.t))
ex.2.t<-ex.2.t[index,]
index<-intersect(Genes$Feature,rownames(ex.3.t))
ex.3.t<-ex.3.t[index,]

##########################
e<-ex.1.t
name<-"1.tumor"
set="GSE14520"
cluster3<-CLUSTER3(e,name)
cluster2<-CLUSTER2(e,name)
cluster3.w<-CLUSTER3.wt(e)
cluster2.w<-CLUSTER2.wt(e)

macro<-data.frame(row.names = colnames(ex.3.t),
                  MFP_M1=as.numeric(e3.mfp[colnames(ex.3.t),"M1_signatures"]),
                  MFP_M=as.numeric(e3.mfp[colnames(ex.3.t),"Macrophages"]),
                  xCell_Monocytes=as.numeric(e3_xcell["Monocytes",colnames(ex.3.t)]),
                  xCell_M0=as.numeric(e3_xcell["Macrophages",colnames(ex.3.t)]),
                  xCell_M1=as.numeric(e3_xcell["Macrophages M1",colnames(ex.3.t)]),
                  xCell_M2=as.numeric(e3_xcell["Macrophages M2",colnames(ex.3.t)]),
                  CIBERSORT_Monocytes=as.numeric(e3_CIBERSORT["Monocytes",colnames(ex.3.t)]),
                  CIBERSORT_M0=as.numeric(e3_CIBERSORT["Macrophages M0",colnames(ex.3.t)]),
                  CIBERSORT_M1=as.numeric(e3_CIBERSORT["Macrophages M1",colnames(ex.3.t)]),
                  CIBERSORT_M2=as.numeric(e3_CIBERSORT["Macrophages M2",colnames(ex.3.t)]),
                  MCP_Monocytic=as.numeric(e3_MCP["Monocytic lineage",colnames(ex.3.t)])
)
cd8t<-data.frame(row.names = colnames(ex.3.t),
                 MFP_T_cell_traffic=as.numeric(e3.mfp[colnames(ex.3.t),"T_cell_traffic"]),
                 MFP_T_cells=as.numeric(e3.mfp[colnames(ex.3.t),"T_cells"]),
                 xCell_CD8_T=as.numeric(e3_xcell["CD8+ T-cells",colnames(ex.3.t)]),
                 xCell_CD8_Tem=as.numeric(e3_xcell["CD8+ Tem",colnames(ex.3.t)]),
                 xCell_CD8_Tn=as.numeric(e3_xcell["CD8+ naive T-cells",colnames(ex.3.t)]),
                 xCell_CD8_Tcm=as.numeric(e3_xcell["CD8+ Tcm",colnames(ex.3.t)]),
                 CIBERSORT_CD8_T=as.numeric(e3_CIBERSORT["T cells CD8",colnames(ex.3.t)]),
                 MCP_CD8_T=as.numeric(e3_MCP["CD8 T cells",colnames(ex.3.t)]))


pcdim<-40
nc<-2
seed=1

ex<-ex.1.t#[rownames(ex.3.t),]#ex.3.t#all[rownames(ex.3.t),]1#exprMatrix.2#
ex["CD74",]->cd74
ex["SPP1",]->spp1
ex["MIF",]->mif

arrydata<-rbind(MIF=mif,SPP1=spp1,CD74=cd74)%>%t()%>%as.data.frame()
arrydata<-log2(arrydata+1)%>%na.omit()%>%as.data.frame()
cor(arrydata[,1:2])

OPT<-function(x,data){
  ndata<-c(0,0,0)%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    s<-mean(data[i:(i+x),1])
    m<-mean(data[i:(i+x),2])
    cd<-mean(data[i:(i+x),3])
    col<-c(s,m,cd)%>%as.data.frame()
    ndata<-cbind(ndata,col)
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  #alist<-list()
  cdt<-c(0,0,0)%>%as.data.frame()
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-data.frame(c(cmat[1,2:3],cmat[2,3]))
    colnames(cc)<-i
    cdt<-cbind(cdt,cc)
    #alist=c(alist,list(ndt.i))
  }
  cdt<-cdt[,-1]%>%t()%>%as.data.frame()
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(4,arrydata)

arrydata<-cbind(t(ex[c("SPP1","MIF"),]),macro,cd8t)%>%as.data.frame()

OPT<-function(x,data){
  ndata<-c(rep(0,ncol(data)))%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  for(i in 1:(nrow(data)-x)){
    col<-c()
    for(j in 1:ncol(data)){
      s<-mean(data[i:(i+x),j])
      col<-c(col,s)
    }
    ndata<-cbind(ndata,as.data.frame(col))
  }
  ndata<-ndata[,-1]
  rownames(ndata)<-colnames(data)
  ndata<-t(ndata)%>%as.data.frame()
  
  return(ndata)
}
OPTy<-function(y,data){
  cdt<-c(rep(0,ncol(data)*2-4))%>%as.data.frame()%>%t()
  colnames(cdt)<-c(paste("SPP1",colnames(data)[-2:-1],sep = "-"),paste("MIF",colnames(data)[-2:-1],sep = "-"))
  for(i in 1:y){
    ndt.i=OPT(i,data)
    cmat<-cor(ndt.i)
    cc<-t(data.frame(c(cmat[1,3:ncol(data)],cmat[2,3:ncol(data)])))
    colnames(cc)<-colnames(cdt)
    cdt<-rbind(cdt,cc)
  }
  cdt<-cdt[-1,]
  rownames(cdt)<-c(2:(y+1))
  return(cdt)
}
arrycdt<-OPTy(30,arrydata)
arrydata.op<-OPT(9,arrydata)

dir.create(name)

if(correlation){
  library(gglmannotate)
  jpeg(file=paste0(name,"\\array-",set,"-MIF-SPP1.jpeg"),width=4,height=4,units = "in", res = 1000)
  ggplot(data=arrydata.op, aes(x=MIF, y=SPP1)) +
    geom_point(alpha=1, size=1,color="#5BBCD6") +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("MIF (logTPM)") + ylab("SPP1 (logTPM)") + 
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
    geom_smooth(aes(x=MIF, y=SPP1),method = glm,linetype=1,se=T,size=1, formula=y ~ x,colour="black",
                fullrange=T)+geom_lmannotate(region = c(xmin = 0, xmax = 0.9, ymin = 0.5, ymax = 1),
                                             place = "topleft")+
    annotate("text",x=10.1,y=6,label="R:0.53",parse=T,size=4)
  dev.off()
  
  bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
  jpeg(file=paste0(name,"\\cor.",set,".infi.jpeg"),width=8,height=6,units = "in", res = 1000)
  pheatmap::pheatmap(t(na.omit(t(arrycdt))),#scale = "row",
                     #kmeans_k = 3,
                     show_colnames =T,
                     show_rownames = T,
                     cluster_row = F, cluster_col = T,
                     border_color ='white',
                     angle_col=45,
                     number_color="gray",
                     display_numbers =round(t(na.omit(t(arrycdt))), 2),
                     fontsize=5,
                     #fontsize_col=10,fontsize_row=5,
                     color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                               colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                     legend_breaks=seq(-1,1,0.2),
                     breaks=bk,
                     main ="correlation optimization") 
  dev.off()
}

if(circle){
  cora<-rbind(SPP1=arrycdt[9,1:(ncol(arrycdt)/2)],MIF=arrycdt[9,(ncol(arrycdt)/2+1):ncol(arrycdt)])%>%as.data.frame()
  cora<-cora[,-7]
  colnames(cora)<-colnames(arrydata)[c(-2:-1,-7)]
  h<-pheatmap::pheatmap(t(cora),scale = "none",
                        show_colnames =T,
                        show_rownames = T,
                        cluster_row = T, cluster_col = F,
                        border_color =NA,
                        number_color="black",
                        display_numbers =round(t(cora), 2),
                        cellwidth = 30, cellheight = 30,
                        angle_col=45,fontsize=15,fontsize_col=12,fontsize_row=12,
                        na_col = "grey",
                        color = c(colorRampPalette(colors = c("#5BBCD6","white"))(length(bk)/2),
                                  colorRampPalette(colors = c("white","#B40F20"))(length(bk)/2)),
                        legend_breaks=seq(-1,1,0.5),
                        breaks=bk,
                        main ="correlation") 
  jpeg(file=paste0(name,"\\cor.",set,".jpeg"),width=6,height=8,units = "in", res = 1000)
  h
  dev.off()
  
  fdf<-t(cora)
  name<-rownames(fdf)[h$tree_row$order]
  fdf<-fdf[name,]
  
  library(reshape3)
  cirmat<-melt(fdf)
  head(cirmat)
  
  resn <- cirmat %>% filter(Var2 == 'MIF')
  rese <- cirmat %>% filter(Var2 == 'SPP1')
  
  resn$ang <- seq(from = (360/nrow(resn)) / 1.5,
                  to = (1.5* (360/nrow(resn))) - 360,
                  length.out = nrow(resn)) + 80
  
  resn$hjust <- 0
  resn$hjust[which(resn$ang < -90)] <- 1
  resn$ang[which(resn$ang < -90)] <- (180+resn$ang)[which(resn$ang < -90)]
  rese$hjust<-resn$hjust
  rese$ang<-resn$ang
  cirmat$var <- rep(1:nrow(resn),2)
  
  s_anno <- data.frame(value = c("SPP1","MIF"))
  s_anno$value <-factor(s_anno$value,levels = c("SPP1","MIF"))
  cirmat$Var2<-factor(cirmat$Var2,levels = c("SPP1","MIF"))
  
  library(gplots)
  library(ggdendro)
  library(ggh4x)
  library(ggnewscale)
  library(ggrepel)
  p1 <- ggplot() +
    geom_bar(data = s_anno,stat = 'identity',
             aes(x = 0,y = 1,fill = value),
             width = 1,
             color = NA) +
    scale_fill_manual(name = 'Gene',
                      values = c("MIF"="#e3D200", 
                                 "SPP1"="#F98400"))
  p2<-p1+ new_scale("fill")+
    geom_tile(data = cirmat[which(cirmat$Var2 == 'MIF'),],
              aes(x = 1:nrow(resn),y = 0.5,fill = value),
              color = 'white') +
    geom_tile(data = cirmat[which(cirmat$Var2 == 'SPP1'),],
              aes(x = 1:nrow(resn),y = 1.5,fill = value),
              color = 'white')+  scale_fill_gradient2(midpoint = 0,
                                                      low = "#5BBCD6",
                                                      mid = "white",
                                                      high = "#B40F20") +
    ylim(-1,5)+xlim(-1,20)
  
  name<-"3.tumor"
  jpeg(file=paste0(name,"\\cor.",set,".circle.jpeg"),width=15,height=10,units = "in", res = 1000)
  p2 + coord_polar(theta = 'x') +theme_void() +new_scale("colour")+
    geom_text(data = resn,
              aes(x = as.numeric(rownames(resn)),
                  y = 2.2,
                  label = as.character(Var1), angle = ang, hjust = hjust),
              size = 4)+
    geom_text(data = resn,
              aes(x = as.numeric(rownames(resn)),
                  y = c(0.25),
                  label = round(value,2), angle = ang, hjust = hjust),
              size = 3.5)+
    geom_text(data = rese,
              aes(x = as.numeric(rownames(rese)),
                  y = c(1.25),
                  label = round(value,2),angle = ang, hjust = hjust),
              size = 4.5)#+NoLegend()
  # annotate("text", x = as.numeric(rownames(rese)), 
  #          y = 0.25, label = round(rese$value,2), color = "white", size = 4)
  dev.off()
}

###############################
library(survival)
library(survivalROC)
library(beeswarm)
library("survminer")
library(preprocessCore)

time<-time.1
#time<-time.2[,c(1,3,2)]
colnames(time)<-c("id","futime", "fustate")

# time<-cbind(time,t(ex.1.t[c("SPP1","MIF"),]))
# time<-cbind(time,t(ex.2.t[c("SPP1","MIF"),]))
# time<-cbind(time.3,t(ex.3.t[c("SPP1","MIF"),]))

time$cluster<-cluster3.w[rownames(time)]
time$fustate<-ifelse(time$fustate=="Alive",0,1)
time$futime<-as.numeric(time$futime)

rt<-time
#rt<-rt[rt$cluster!="CT",]
rt<-na.omit(rt)
rt$futime=rt$futime#*365/30
table(rt$cluster)
max(rt$futime)

#rt<-rt[rt$futime<=60,]

fit <- survfit(Surv(futime, fustate) ~ cluster, data = rt)
diff=survdiff(Surv(futime, fustate) ~cluster,data = rt)
pValue=1-pchisq(diff$chisq,df=length(diff$n)-1)
pValue=signif(pValue,4)
pValue
#0.0498
#0.08209
colsurv=c("CT"="#FF0000","AN"="#00A08A","PT"="#F98400")
jpeg(file=paste0(name,"\\0907-survival-",set,"-intergrated.jpeg"),width =5,height =6,units = "in", res = 2000)
ggsurvplot(fit, 
           data=rt,
           #conf.int=TRUE,
           palette=c("#00A08A","#F98400","#FF0000"),#c("#F98400","#E2D200","#00A08A","#5BBCD6"),#
           pval=paste0("p=",pValue),
           pval.size=4,
           #legend.labs=c("High", "Low"),
           #legend.title="Expression",
           break.time.by = 6,
           risk.table.title="",
           risk.table=TRUE,
           risk.table.height=.30,
           #xlim=c(0,24),
           ylim=c(0,1),
           xlab="Time(months)")
dev.off()

save(df,time.1,ex.1.t,file = "D:\\Bioinfrolf\\Bioinfrolf\\tmpRdata/GSE14520_intergrated.Rdata")

90+67+50
df1<-na.omit(time)
max(df1$futime)

#df1$futime=df1$futime/30
df1<-df1[order(df1$MIF,decreasing = T),]
df1$levels<-ifelse(df1$MIF>=median(df1$MIF),"High expression","Low expression")

df1<-df1[df1$futime<=60,]
outTab=data.frame()
for (i in 10:60){
  for(j in 10:60){
    rt<-df1[c(1:i,(nrow(df1)-j+1):nrow(df1)),]
    fit <- survfit(Surv(futime, fustate) ~ levels, data = rt)
    diff=survdiff(Surv(futime, fustate) ~levels,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    pValue=signif(pValue,4)
    out=data.frame(h=i,l=j,p=pValue)
    outTab=rbind(outTab,out)
  }
}

outTab<-outTab[order(outTab$p,decreasing = F),]
outTab$d<-abs(outTab$h-outTab$l)
outTab$a<-outTab$h+outTab$l
outTab<-outTab[outTab$p<0.05,]
outTab.1<-outTab[outTab$d<50,]
outTab.1<-outTab.1[ outTab.1$a <171,]
outTab.1<-outTab.1[order(outTab.1$p,decreasing = F),]
write.csv(outTab,"MIF-",set,"pvalue.csv")

i=84
j=82
rt<-df1[c(1:i,(nrow(df1)-j+1):nrow(df1)),]
fit <- survfit(Surv(futime, fustate) ~ levels, data = rt)
diff=survdiff(Surv(futime, fustate) ~levels,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)

jpeg(file=paste0(name,"\\survival-",set,"-5years-MIF-1.jpeg"),width =5,height =6,units = "in", res = 1000)
ggsurvplot(fit, 
           data=rt,
           #conf.int=TRUE,
           palette=c("#B40F20","#5BBCD6"),
           pval=paste0("p=",pValue),
           pval.size=4,
           legend.labs=c("High", "Low"),
           legend.title="Expression",
           break.time.by = 6,
           risk.table.title="",
           risk.table=TRUE,
           risk.table.height=.30,
           #xlim=c(0,24),
           ylim=c(0.2,1),
           xlab="Time(months)")
dev.off()

save(list = ls(),file = "TCGA.Rdata")