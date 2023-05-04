library(limma)
library(umap)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

rm(list = ls())
setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\prognosis")
stDEs=read.table("ME&MI-DE-seq.txt",header = T,row.names = 1)
sDEs=read.table("T&N-DE-seq.txt",header = T,row.names = 1)

setwd("D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge")
stDEa=read.table("ME&MI-DE-1.txt",header = T,row.names = 1)
sDEa=read.table("T&N-DE-array.txt",header = T,row.names = 1)

stDEs$p=ifelse(stDEs$pValue<0.05,1,0)
stDEa$p=ifelse(stDEa$pValue<0.05,1,0)
sDEs$p=ifelse(sDEs$pValue<0.05,1,0)
sDEa$p=ifelse(sDEa$pValue<0.05,1,0)

un=union(rownames(stDEs),rownames(stDEa))
dimnames=list(un,c("seq","array"))
undf=matrix(rep(0,length(un)*2),nrow=length(un),dimnames=dimnames)%>%as.data.frame()
undf[rownames(stDEs),1]=stDEs[,"logFC"]
undf[rownames(stDEa),2]=stDEa[,"logFC"]

inm=intersect(rownames(stDEs),rownames(stDEa))
dimnames=list(inm,c("seq","array"))
inmdf=matrix(rep(0,length(inm)*2),nrow=length(inm),dimnames=dimnames)%>%as.data.frame()
inmdf$seq=stDEs[inm,"logFC"]
inmdf$array=stDEa[inm,"logFC"]

intersect(rownames(stDEs)[stDEs$p==1],rownames(stDEa)[stDEa$p==1])
inmdf[which(rownames(inmdf)%in%c("MIF","CD74","SPP1")),]

inmdf$result = as.factor(ifelse(inmdf$seq>0.5 & inmdf$array > 0.1,"high risk",
                                ifelse(inmdf$seq<(-0.5) & inmdf$array<(-0.1),
                                       "low risk",'not')))
inmdf$Label = c(rep(NA,nrow(inmdf)))
inmdf$Label[which(rownames(inmdf)%in%c("MIF","SPP1"))]=c("MIF","SPP1")
co=cor(as.numeric(inmdf$seq),as.numeric(inmdf$array))

int=intersect(rownames(sDEs),rownames(sDEa))
dimnames=list(int,c("seq","array"))
intdf=matrix(rep(0,length(int)*2),nrow=length(int),dimnames=dimnames)%>%as.data.frame()
intdf$seq=sDEs[int,"logFC"]
intdf$array=sDEa[int,"logFC"]

intersect(rownames(sDEs)[sDEs$p==1],rownames(sDEa)[sDEa$p==1])
intdf[which(rownames(intdf)%in%c("MIF","CD74","SPP1")),]

intdf$result = as.factor(ifelse(intdf$seq>0.5 & intdf$array > 0.1,"high risk",
                                ifelse(intdf$seq<(-0.5) & intdf$array<(-0.1),
                                       "low risk",'not')))
intdf$Label = c(rep(NA,nrow(intdf)))
intdf$Label[which(rownames(intdf)%in%c("MIF","SPP1"))]=c("SPP1","MIF")
co=cor(as.numeric(intdf$seq),as.numeric(intdf$array))

hr<-list(t=rownames(intdf)[intdf$result=="high risk"],
        m=rownames(inmdf)[inmdf$result=="high risk"])
intersect(hr$t,hr$m)
lr<-list(t=rownames(intdf)[intdf$result=="low risk"],
         m=rownames(inmdf)[inmdf$result=="low risk"])
intersect(lr$t,lr$m)

library(VennDiagram)   
library(gplots)
col2<-c("#B40F20","#B40F20")
col2=c("#5BBCD6","#5BBCD6")
barplot(1:8,col=col8)
?venn.diagram()
venn.diagram(
  x = hr,
  filename = 'VN-h.png',
  col = "black",
  fill = col2,
  alpha = 0.4,
  cex = 4,
  cat.col = 'black',
  cat.cex = 0,
  cat.fontface = "bold",
  margin = 0.1,
  main = "high risk genes",
  main.cex = 5
)

jpeg(file="Tumor-DE.jpeg",width=6,height=4,units = "in", res = 600)
ggplot(intdf,aes(x=array,y=seq))+geom_point(aes(color = result),size=2)+theme_bw()+
  theme(panel.grid = element_blank(),
                   axis.line = element_line(size = 0.5),
                   axis.text.x =element_text(face = 'bold',size = 10),
                   axis.ticks.x = element_blank(),
                   axis.text.y = element_text(face = 'bold',size = 10),
  )+xlab('logFC(array)') + ylab('logFC(seq)')+ggtitle("DE")+
  scale_color_manual(values =c("#B40F20","#5BBCD6","gray"))+
  geom_hline(yintercept=c(-0.5,0.5),linetype="dashed",size=0.5)+ 
  geom_vline(xintercept= c(-0.1,0.1),linetype="dashed",size=0.5)+
  geom_label_repel(aes(label=Label), fontface="bold",size=3,
                   box.padding=unit(0.2, "lines"), point.padding=unit(0.3, "lines"),
                   segment.colour = "black",na.rm = TRUE,show.legend = FALSE)
dev.off()
#df=data.frame(seqTN=sDEs,seqM=stDEs,arrayTN=sDEa,arrayM=stDEa)



write.csv(intdf,"limma_results-20210603.csv",quote = F)
###########vol-1#################
this_tile2 <- paste0('Cutoff for logFC is ',round(logFC_cutof,3), #round????С??λ??
                     '\nThe number of up gene is ',nrow(intdf[intdf$result =='UP',]) ,
                     '\nThe number of down gene is ',nrow(intdf[intdf$result =='DOWN',]))

p2<-ggplot(data=intdf, aes(x=logFC, y=-log10(P.Value), color=result)) +
  geom_point(alpha=0.4, size=2) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value")+
  ggtitle( this_tile2 ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c("darkgreen", "grey", "brown4"),limits = c("DOWN","NOT","UP"))+
  theme_bw()+theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2<-p2+ geom_hline(yintercept=-log10(fdrFilter),linetype="dashed",size=1)+ 
  geom_vline(xintercept= c(-logFC_cutof,logFC_cutof),linetype="dashed",size=1)+
  theme(axis.text.x=element_text(hjust =0.5, vjust =0.5,size=18,face = "bold"),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=20,face = "bold"),
        axis.title.x = element_text(size=20,face = "bold"))
#intdf$Label <- ifelse(intdf$adj.P.Val < 0.001 & abs(intdf$logFC) > 2.5,rownames(intdf),NA)
p3<-p2+geom_label_repel(aes(label=intdf$Label), 
                        fontface="bold",
                        box.padding=unit(0.2, "lines"), 
                        point.padding=unit(0.3, "lines"), 
                        segment.colour = "grey50",
                        na.rm = TRUE)