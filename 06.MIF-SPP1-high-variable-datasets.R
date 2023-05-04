library(dplyr)
library(ggplot2)
load("/Users/zhaolab3/Documents/scRNA_MNP/oldAnalysis/ICGC_count.Rdata")
sample<-meta$submitted_sample_id
tissue<-unlist(lapply(sample,function(x){
  unlist(strsplit(x,"_"))[2]
}))
tissue<-ifelse(tissue=="Liver","Normal","Tumor")
table(tissue)
meta$tissue<-tissue
rm(data.raw)

setwd("/Users/zhaolab3/Documents/scRNA_MNP/")
tissue=c(rep(c("Normal","Tumor"),c(50,374)),meta$tissue)
batch=rep(c("TCGA","ICGA"),c(424,445))

fdrFilter<-0.05
conNum=50           
treatNum=374
exprMatrix.2 = read.table(file = "./oldAnalysis/TCGA_symbol_tpm.txt",
                          header=TRUE,row.names=1, as.is=TRUE)

exprMatrix.3 = read.table(file = "./oldAnalysis/ICGC_symbol_tpm.txt",
                          header=TRUE,row.names=1, as.is=TRUE)
exprMatrix.3[1:5,1:5]

ex.TCGA.t<-exprMatrix.2[,51:ncol(exprMatrix.2)]
ex.ICGC.t<-exprMatrix.3[,meta$tissue=="Tumor"]

ex.TCGA.t[1:5,1:5]

l<-intersect(rownames(ex.TCGA.t),rownames(ex.ICGC.t))
###################################
OD<-function(data,gene){
  data<-t(data)%>%as.data.frame()
  index<-which(colnames(data)==gene)
  data<-data[order(data[,index],decreasing = F),]
  data<-t(data)%>%as.data.frame()
}
gene<-"MIF"
ex.TCGA.t<-OD(ex.TCGA.t,gene)
ex.ICGC.t<-OD(ex.ICGC.t,gene)

outTab=data.frame()
name<-"TCGA-"
rt<-ex.TCGA.t
j=50
for (j in seq(100,150,50)){
  treatNum=j
  conNum=j 
  outTab=data.frame()
  grade=c(rep(1,conNum),rep(2,treatNum))
    data<-rt[,c(1:j,(ncol(rt)-j+1):ncol(rt))]
    for(i in row.names(data)){
      tmp=rbind(expression=data[i,],grade=grade)
      tmp=as.matrix(t(tmp))
      wilcoxTest<-wilcox.test(expression ~ grade, data=tmp)
      conGeneMeans=mean(as.numeric(data[i,1:conNum]))
      treatGeneMeans=mean(as.numeric(data[i,(conNum+1):ncol(data)]))
      logFC=log2(treatGeneMeans)-log2(conGeneMeans)
      pvalue=wilcoxTest$p.value
      #conMed=median(as.numeric(data[i,1:conNum]))
      #treatMed=median(as.numeric(data[i,(conNum+1):ncol(data)]))
      #diffMed=treatMed-conMed
      #logFC.med=log2(treatMed)-log2(conMed)
      #if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){   } 
      outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,
                                #conMed=conMed,treatMed=treatMed,logFC.med=logFC.med,
                                pValue=pvalue))
      
    }
    pValue=outTab[,"pValue"]
    fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
    outTab=cbind(outTab,fdr=fdr)
    outTab$pValue<-as.numeric(outTab$pValue)
    outTab=outTab[order(outTab$pValue,decreasing = F),]
    res<-outTab[which(outTab$gene%in%c("MIF","CD74","SPP1")),]
    write.csv(outTab,file = paste0(name,gene,"-",j,"-DEG.csv"),quote = F)
    write.csv(res,file = paste0(name,gene,"-",j,"-res.csv"),quote = F)
}
#################DEG#################
j=50
IM<-read.csv("./oldAnalysis/ICGC-MIF-50-DEG.csv")
IS<-read.csv("./oldAnalysis/ICGC-SPP1-50-DEG.csv")
TM<-read.csv("./oldAnalysis/TCGA-MIF-50-DEG.csv")
TS<-read.csv("./oldAnalysis/TCGA-SPP1-50-DEG.csv")
G2<-c("MIF","SPP1")

gene<-"MIF"#
name<-"ICGC-"#"TCGA-"
which(IM$gene=="SPP1")

nrDEG<-IM
nrDEG$logFC<-as.numeric(nrDEG$logFC)
nrDEG$logFC[which(nrDEG$logFC%in%c(Inf,-Inf))]=0
nrDEG$pValue<-as.numeric(nrDEG$pValue)
nrDEG<-nrDEG[order(nrDEG$pValue,decreasing = F),]
nrDEG<-na.omit(nrDEG)
logFC_cutof = with(nrDEG,mean(abs(nrDEG$logFC)) + 2*sd(abs(nrDEG$logFC)))
logFC_cutof 
#logFC_cutof=nrDEG$logFC[which(nrDEG$gene=="MIF")]
nrDEG[which(nrDEG$gene%in%G2),]
nrDEG$result = as.factor(ifelse(nrDEG$pValue < 0.05 & abs(nrDEG$logFC) >=logFC_cutof,
                                ifelse(nrDEG$logFC >= logFC_cutof ,'UP','DOWN'),'NOT'))
# nrow(nrDEG[nrDEG$result =='UP',]) 
# nrow(nrDEG[nrDEG$result =='DOWN',])
# up.genes<- head(nrDEG$gene[which(nrDEG$result =="UP")],10)
# down.genes<- head(nrDEG$gene[which(nrDEG$result =="DOWN")],10)
# deg.top10.genes<- c(as.character(up.genes),as.character(down.genes))  
rownames(nrDEG)<-nrDEG$gene
nrDEG$Label = c(rep(NA,nrow(nrDEG)))
fc<-signif(nrDEG$logFC[match(G2,rownames(nrDEG))],3)
nrDEG$Label[match(G2,rownames(nrDEG))]<- paste(G2,fc,sep = " LogFC:")

this_tile2 <- paste0('Cutoff for logFC is ',round(logFC_cutof,3))

p2<-ggplot(data=nrDEG, aes(x=logFC, y=-log10(pValue), color=result)) +
  geom_point(alpha=0.4, size=0.5) +
  theme_set(theme_set(theme_bw(base_size=10)))+
  xlab("log2 fold change") + ylab("-log10 p-value")+
  ggtitle( this_tile2 ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c("#5BBCD6","grey","#B40F20"),limits = c("DOWN","NOT","UP"))+
  theme_bw()+theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2<-p2+ geom_hline(yintercept=-log10(fdrFilter),linetype="dashed",size=0.5)+ 
  geom_vline(xintercept= c(-logFC_cutof,logFC_cutof),linetype="dashed",size=0.5)+
  theme(axis.text.x=element_text(hjust =0.5, vjust =0.5,size=12),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12))

jpeg(paste0("vol-",name,gene,"-",j,"-DEG.jpg"),width =4.5,height = 3,units = "in", res = 1000)
p2+geom_label_repel(aes(label=nrDEG$Label), 
                    fontface="bold",
                    #box.padding = F,
                    #box.padding=unit(0.2, "lines"), 
                    #point.padding=unit(0.3, "lines"), 
                    #segment.colour = "grey50",
                    na.rm = TRUE)
dev.off()
#################Venn#################################
G2<-c("MIF","SPP1")

gene<-"MIF"#
name<-"ICGC-"#"TCGA-"
which(IM$gene=="SPP1")
filterRES<-function(res){
  res<-na.omit(res)
  res<-res[abs(res$logFC)!=Inf,]
  res<-res[res$pValue<0.05,]
  res<-res[order(res$logFC,decreasing = T),]
  rownames(res)<-res$gene
  return(res)
}
IM<-filterRES(IM)
IS<-filterRES(IS)
TM<-filterRES(TM)
TS<-filterRES(TS)
MU<-intersect(IM$gene[IM$logFC>0],TM$gene[TM$logFC>0])
MD<-intersect(IM$gene[IM$logFC<0],TM$gene[TM$logFC<0])
SU<-intersect(IS$gene[IS$logFC>0],TS$gene[TS$logFC>0])
SD<-intersect(IS$gene[IS$logFC<0],TS$gene[TS$logFC<0])

hr<-list(g1=MU,g2=SU)
h<-intersect(hr$g1,hr$g2)
Muni<-MU[!MU%in%h]
Suni<-SU[!SU%in%h]
lr<-list(g1=MD,g2=SD)
l<-intersect(lr$g1,lr$g2)

library(VennDiagram)   
library(gplots)
col3<-c("#E2D200","#00A08A")

col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6",
        "#E2D200","#B40F20","#273046","#FD6467")
barplot(1:8,col=col8)
library(eulerr)
vd <- euler(c("MIF Up"= (length(hr$g1)-3049),
              "SPP1 Up"=(length(hr$g2)-3049),
              "MIF Up&SPP1 Up" = 3049))
pvalue <- phyper(3049-1,7513,17241-7513,4197,lower.tail=F)
qvalue <- p.adjust(pvalue,method='fdr')
qvalue

jpeg(file="MU_SU1.jpeg",width=4,height=1,units = "in", res = 1000)
plot(vd,
     fills = list(fill =c("#B40F20","#FF0000"), alpha = 0.4),
     labels = list(col = "black", font = 1,size=0.1), 
     edges = T,
     quantities = TRUE)
dev.off()

vd <- euler(c("MIF Down"= (length(lr$g1)-442),
             "SPP1 Down"=(length(lr$g2)-442),
              "MIF Down&SPP1 Down" = 442))
phyper(442-1,799,17241-799,784,lower.tail = F)%>%p.adjust()

jpeg(file="MD_SD1.jpeg",width=4,height=1,units = "in", res = 1000)
plot(vd,
     fills = list(fill = c("#5BBCD6","#00A08A"), alpha = 0.4),
     labels = list(col = "black", font = 1), 
     edges = T,
     quantities = TRUE)
dev.off()

save(list = ls(),file = "tmp.Rdata")
library(presto)
head(pbmc.genes)
#################################
library(tibble)
library(msigdbr)
library(fgsea)
hs_msigdbr <- msigdbr(species="Homo sapiens")
hs_msigdbr %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
hsKEGG <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
hsGO <- msigdbr(species="Homo sapiens",category="C5")
#hsGOBP <- msigdbr(species="Homo sapiens",category="C5",subcategory = "BP")

fgsea_sets<- hsGO %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_KEGG<-hsKEGG %>% split(x = .$gene_symbol, f = .$gs_name)
cluster<-rep(c("high","low"),c(50,50))

runGSEA<-function(nrDEG){
  cluster.genes<- nrDEG%>% arrange(desc(logFC)) %>% dplyr::select(gene, logFC)
  ranks<- deframe(cluster.genes)
  library(msigdbr)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
  fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
  return(fgseaResTidy)
}
runGSEA_kegg<-function(nrDEG){
  cluster.genes<- nrDEG%>% arrange(desc(logFC)) %>% dplyr::select(gene, logFC)
  ranks<- deframe(cluster.genes)
  library(msigdbr)
  fgseaRes<- fgsea(fgsea_KEGG, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
  fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
  return(fgseaResTidy)
}
keggIM<-runGSEA_kegg(IM)
goIM<-runGSEA(IM)
keggIS<-runGSEA_kegg(IS)
goIS<-runGSEA(IS)
keggTS<-runGSEA_kegg(TS)
goTS<-runGSEA(TS)
keggTM<-runGSEA_kegg(TM)
goTM<-runGSEA(IM)
filterRES<-function(res){
  res<-na.omit(res)
  res<-res[res$pval<0.05,]
  res<-res[order(res$NES,decreasing = T),]
  res<-as.data.frame(res)
  rownames(res)<-res$pathway
  res$pathway<-lapply(res$pathway,function(x){
    paste(unlist(strsplit(x,split = "_"))[-1],
          collapse = " ")%>% str_to_title()
  })%>%unlist()
  a<-gsub("In","in",res$pathway)
  a<-gsub("Of","of",a)
  res$pathway<-gsub("To","to",a)
  return(res)
}
goIM<-filterRES(goIM)
goTM<-filterRES(goTM)
goIS<-filterRES(goIS)
goTS<-filterRES(goTS)
copMg<-intersect(rownames(goTM),rownames(goIM))
copSg<-intersect(rownames(goTS),rownames(goIS))

keggIM<-filterRES(keggIM)
keggTM<-filterRES(keggTM)
keggIS<-filterRES(keggIS)
keggTS<-filterRES(keggTS)

copMk<-intersect(rownames(keggTM),rownames(keggIM))
copSk<-intersect(rownames(keggTS),rownames(keggIS))
keggIM$set<-"ICGC"
keggTM$set<-"TCGA"
keggIS$set<-"ICGC"
keggTS$set<-"TCGA"

goIM$set<-"ICGC"
goTM$set<-"TCGA"
goIS$set<-"ICGC"
goTS$set<-"TCGA"
MG<-rbind(goIM[copMg,c(1,2,5,9)],goTM[copMg,c(1,2,5,9)])
SG<-rbind(goIS[copSg,c(1,2,5,9)],goTS[copSg,c(1,2,5,9)])

{
  MK<-rbind(keggIM[copMk,c(1,2,5,9)],keggTM[copMk,c(1,2,5,9)])
  mk<-cbind(keggIM[copMk,c(1,2,5,9)],keggTM[copMk,c(1,2,5,9)])
  mk$NESA<-mg[,3]+mg[,7]
  MK$posation<-MK$NES
  MK$posation[1:27]<-mk$NESA
  MK<-MK[order(MK$posation,decreasing = T),]
  MK$pathway<-factor(MK$pathway,levels = unique(MK$pathway))
  MK<-MK[MK$pathway!="Maturity Onset Diabetes of The Young",]
  jpeg(file="MIF_KEGG-1.jpeg",width =7,height =5,units = "in", res = 1000)
  ggplot(MK,aes(x=NES,y=pathway,fill=set))+geom_bar(stat="identity",position="stack",width=0.2)+
    geom_point(aes(x=posation,color =pval), size = 4)+
    scale_color_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
    scale_fill_manual(values=c("#5BBCD6","#E2D200"))+
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
  
  SK<-rbind(keggIS[copSk,c(1,2,5,9)],keggTS[copSk,c(1,2,5,9)])
  sk<-cbind(keggIS[copSk,c(1,2,5,9)],keggTS[copSk,c(1,2,5,9)])
  sk$NESA<-sk[,3]+sk[,7]
  SK$posation<-SK$NES
  SK$posation[1:30]<-sk$NESA
  SK<-SK[order(SK$posation,decreasing = T),]
  SK$pathway<-factor(SK$pathway,levels = unique(SK$pathway))
  jpeg(file="SPP1_KEGG-1.jpeg",width =6.65,height =5.5,units = "in", res = 1000)
  ggplot(SK,aes(x=NES,y=pathway,fill=set))+geom_bar(stat="identity",position="stack",width=0.2)+
    geom_point(aes(x=posation,color =pval), size = 4)+
    scale_color_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
    scale_fill_manual(values=c("#5BBCD6","#E2D200"))+
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
}

MG<-rbind(goIM[copMg,c(1,2,5,9)],goTM[copMg,c(1,2,5,9)])
mg<-cbind(goIM[copMg,c(1,2,5,9)],goTM[copMg,c(1,2,5,9)])
mg$NESA<-mg[,3]+mg[,7]
mg<-mg[abs(mg$NESA)>mg[,3]&&abs(mg$NESA)>mg[,7],]
MG$posation<-MG$NES
MG$posation[1:1268]<-mg$NESA
MG<-MG[order(MG$posation,decreasing = T),]
MG$pathway<-factor(MG$pathway,levels = unique(MG$pathway))

SG<-rbind(goIS[copSg,c(1,2,5,9)],goTS[copSg,c(1,2,5,9)])
sg<-cbind(goIS[copSg,c(1,2,5,9)],goTS[copSg,c(1,2,5,9)])
sg$NESA<-sg[,3]+sg[,7]
sg<-sg[abs(sg$NESA)>sg[,3]&&abs(sg$NESA)>sg[,7],]
SG$posation<-SG$NES
SG$posation[1:1239]<-sg$NESA
SG<-SG[order(SG$posation,decreasing = T),]
SG$pathway<-factor(SG$pathway,levels = unique(SG$pathway))

cog<-intersect(rownames(mg),rownames(sg))
mg<-mg[cog,]
sg<-sg[cog,]

MG<-rbind(goIM[cog,c(1,2,5,9)],goTM[cog,c(1,2,5,9)])
MG$posation<-MG$NES
MG$posation[1:506]<-mg$NESA
mg<-mg[order(mg$NESA,decreasing = T),]
MG$pathway<-factor(MG$pathway,levels = unique(mg$pathway))
MG<-MG[MG$pathway%in%mg$pathway[c(1:15,492:506)],]

jpeg(file="MIF_go.jpeg",width =10,height =5,units = "in", res = 1000)
ggplot(MG,aes(x=NES,y=pathway,fill=set))+geom_bar(stat="identity",position="stack",width=0.2)+
  geom_point(aes(x=posation,color =pval), size = 4)+
  scale_color_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
  scale_fill_manual(values=c("#5BBCD6","#E2D200"))+
  scale_y_discrete(labels=function(y) str_wrap(y,width = 50))+
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

SG<-rbind(goIS[copSg,c(1,2,5,9)],goTS[copSg,c(1,2,5,9)])
sg<-cbind(goIS[copSg,c(1,2,5,9)],goTS[copSg,c(1,2,5,9)])
sg$NESA<-sg[,3]+sg[,7]
sg<-sg[abs(sg$NESA)>sg[,3]&&abs(sg$NESA)>sg[,7],]
SG$posation<-SG$NES
SG$posation[1:1239]<-sG$NESA
SG<-SG[order(SG$posation,decreasing = T),]
SG$pathway<-factor(SG$pathway,levels = unique(SG$pathway))
jpeg(file="SPP1_go.jpeg",width =10,height =5.5,units = "in", res = 1000)
ggplot(SG,aes(x=NES,y=pathway,fill=set))+geom_bar(stat="identity",position="stack",width=0.2)+
  geom_point(aes(x=posation,color =pval), size = 4)+
  scale_color_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
  scale_fill_manual(values=c("#5BBCD6","#E2D200"))+
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

jpeg(file="SPP1_go.jpeg",width =10,height =5.5,units = "in", res = 1000)
ggplot(SG,aes(x=NES,y=pathway,fill=set))+geom_bar(stat="identity",position="stack",width=0.2)+
  geom_point(aes(x=posation,color =pval), size = 4)+
  scale_color_gradient2(low="#B40F20" ,mid = "white",high="#5BBCD6",midpoint = 0.05)+
  scale_fill_manual(values=c("#5BBCD6","#E2D200"))+
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
cok<-intersect(MK$pathway,SK$pathway)
save(list=ls(),file="MIF_SPP1_GSEA.Rdata")
#####################################################
corFilter=0.4              #????ϵ?????˱?׼
pvalueFilter=0.05         #pֵ???˱?׼

setwd("D:/Bioinfrolf/new figs/")
TF=read.table("files/TF.txt",sep="\t",header=F)
gene<-c(h,l)
outTab=data.frame()
filterFC<-function(nrDEG){
  logFC_cutof = with(nrDEG,mean(abs(nrDEG$logFC)) + 2*sd(abs(nrDEG$logFC)))
  nrDEG$result = as.factor(ifelse(nrDEG$pValue < 0.05 & abs(nrDEG$logFC) >=logFC_cutof,
                                  ifelse(nrDEG$logFC >= logFC_cutof ,'UP','DOWN'),'NOT'))
  g<-nrDEG$gene[nrDEG$result!="NOT"]
  return(g)
}
img<-filterFC(IM)
isg<-filterFC(IS)
tmg<-filterFC(TM)
tsg<-filterFC(TS)

ig<-Reduce(intersect,list(img,isg),accumulate =FALSE)
tg<-Reduce(intersect,list(tmg,tsg),accumulate =FALSE)
cg<-Reduce(intersect,list(ig,tg),accumulate =FALSE)

rt<-ex.TCGA.t
tf=intersect(rownames(rt),TF$V1)

datatf=rt[tf,]
datadg=rt[cg,]

outTab=data.frame()
for(i in row.names(datatf)){
  if(sd(datatf[i,])>1){
    for(j in row.names(datadg)){
      x=as.numeric(datatf[i,])
      y=as.numeric(datadg[j,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,datadg=j,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(TF=i,datadg=j,cor,pvalue,Regulation="negative"))
      }
    }
  }
}
write.table(file="0906_corResult-TCGA-co.txt",outTab,sep="\t",quote=F,row.names=F)        #?????????Խ???


iup<-intersect(IM$gene[IM$logFC>0],IS$gene[IS$logFC>0])
idown<-intersect(IM$gene[IM$logFC<0],IS$gene[IS$logF<0])
tup<-intersect(TM$gene[IM$logFC>0],TS$gene[TS$logFC>0])
tdown<-intersect(TM$gene[IM$logFC<0],TS$gene[TS$logFC<0])

TFLabel=cbind(tf,"TF")

datadgupLabel=cbind(cg[cg%in%tup],"UP")
datadgdownLabel=cbind(cg[cg%in%tdown],"DOWN")
nodeLabel=rbind(c("ID","type"),TFLabel,datadgupLabel,datadgdownLabel)
write.table(nodeLabel,file="nodeType-TCGA-co.txt",sep="\t",quote=F,col.names=F,row.names=F)

datadgupLabel=cbind(cg[cg%in%iup],"UP")
datadgdownLabel=cbind(cg[cg%in%idown],"DOWN")
nodeLabel=rbind(c("ID","type"),TFLabel,datadgupLabel,datadgdownLabel)
write.table(nodeLabel,file="nodeType-ICGC-co.txt",sep="\t",quote=F,col.names=F,row.names=F)

nodeLabel1<-nodeLabel
outTab1<-outTab
outTab1$co<-paste(outTab1$TF,outTab1$datadg,sep = "-")
outTab$co<-paste(outTab$TF,outTab$datadg,sep = "-")
tcor<-outTab1[outTab1$co%in%intersect(outTab1$co,outTab$co),-6]
icor<-outTab[outTab$co%in%intersect(outTab1$co,outTab$co),-6]
write.table(tcor,file="nodeType-TCGA-co-1.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(icor,file="nodeType-ICGC-co-1.txt",sep="\t",quote=F,col.names=T,row.names=F)

########################################
rm(list = ls())
TF=read.table("files/TF.txt",sep="\t",header=F)
Tn<-read.table("nodeType-TCGA-co-1.txt",header = T)

edges <-Tn[,1:2]
colnames(edges) <- c("from", "to")
nodes <- data.frame(name = unique(union(edges$from, edges$to)))e
nodes$type=ifelse(nodes$name%in% TF$V1,"TF","target")
# Tgraph1<-tbl_graph(nodes = nodes, edges = edges)

library(tidygraph)
Tgraph <- as_tbl_graph(Tn, directed = F,edges =Tn$cor)
Igraph <- as_tbl_graph(In, directed = F,edges =In$cor)

centrality_degree(
  weights = NULL,
  mode = "out",
  loops = TRUE,
  normalized = FALSE
)

Tg<-Tgraph %>%
#  activate(Tgraph) %>%
  mutate(centrality = centrality_degree(weights = cor))

jpeg(file="Tn-MIF_KEGG-degree3.jpeg",width =7,height =5,units = "in", res = 1000)
ggraph(Tg,layout = 'linear', circular = TRUE) + 
  geom_edge_link(aes(edge_width=cor),colour="gray") +
  geom_node_point(aes(size = centrality,fill=factor(nodes$type),colour=factor(nodes$type)),shape=21) +
  geom_node_text(aes(filter= centrality>3,label = name),size=4, repel = TRUE)+
  theme_graph()+scale_edge_width(range=c(0.25,1.25))+
  scale_fill_manual(values=c("#F98400","#00A08A"))+
  scale_color_manual(values =c("#F98400","#00A08A"))
dev.off()

Ig<-Igraph %>%
  #  activate(Tgraph) %>%
  mutate(centrality = centrality_degree(weights = cor))

jpeg(file="In-MIF_KEGG-degree3.jpeg",width =7,height =5,units = "in", res = 1000)
ggraph(Ig,layout = 'linear', circular = TRUE) + 
  geom_edge_link(aes(edge_width=cor),colour="gray") +
  geom_node_point(aes(size = centrality,fill=factor(nodes$type),colour=factor(nodes$type)),shape=21) +
  geom_node_text(aes(filter= centrality>3,label = name),size=4, repel = TRUE)+
  theme_graph()+scale_edge_width(range=c(0.25,1.25))+
  scale_fill_manual(values=c("#F98400","#00A08A"))+
  scale_color_manual(values =c("#F98400","#00A08A"))
dev.off()
# factor(Tn$Regulation)
layout = "graphopt"

