#################destinyPseudotime###############
load("D:/Bioinfrolf/Bioinfrolf/destiny/destinyPseudotime.Rdata")
index=intersect(rownames(d2),rownames(d3))
dsdf<-data.frame(row.names = index,Monocle2=d2[index,1],#Monocle2r=d2_r[index,1],
                 Monocle3=d3[index,1],#Monocle3_r=d3_r[index,1],
                 slingshot=ds[index,1],#slingshot_r=ds_r[index,1],
                 destiny=dd[index,1]#,destiny_r=dd_r[index,1]
)
cordf<-cor(dsdf)

index=rownames(d2)
dadf<-data.frame(row.names = index,Monocle2=d2[index,1],#Monocle2r=d2_r[index,1],
                 Monocle3=d3_a[index,1],#Monocle3_r=d3__a_r[index,1],
                 slingshot=ds[index,1],#slingshot_r=ds_r[index,1],
                 destiny=dd[index,1]#,destiny_r=dd_r[index,1]
)
cordfa<-cor(dadf)
library(corrplot)
jpeg(file="corTI-allcell.jpeg",width=7,height=5.5,units = "in", res = 1000)
corrplot(cordfa, type = "upper", tl.pos = "tp",
         tl.cex=1.5,col = colorRampPalette(c("#46ACC8","white","#B40F20"))(50),
         cl.cex = 0.8,#order = "hclust", addrect = 2,
         tl.col = "black")
corrplot(cordfa, add = TRUE, type = "lower", method = "number",
         col = colorRampPalette(c("#46ACC8","white","#B40F20"))(50), 
         cl.cex = 2,
         diag = FALSE, tl.pos = "n", cl.pos = "n")
dev.off()
#################VN########
load("Boruta-noCCL5_2.Rdata")
library(Boruta)
str(feature.selection)
feature.selection%>%str()
jpeg(file="Boruta.jpeg",width=10,height=8,units = "in", res = 1000)
plot(feature.selection)
dev.off()
jpeg(file="Boruta.ImpHistory.jpeg",width=10,height=8,units = "in", res = 1000)
plotImpHistory(feature.selection)
dev.off()

setwd("D:/Bioinfrolf/Bioinfrolf/model")
velo<-read.csv("veloGenes.csv",header = T)
velo<-velo[velo$fit_likelihood>0,]
slingshot<-read.table("slingshot.txt",header = T)
Monocle3<-read.table("Monocle3.gene_noCCL5_all.txt",header = T)
Monocle3<-na.omit(Monocle3)
#macro<-read.table("marcophage genes.txt",header = T)
index<-intersect(Monocle3$gene_short_name,velo$Gene)
Genes<-union(Monocle3$gene_short_name,slingshot$x)

fNamesAd<-names(feature.selection$finalDecision[feature.selection$finalDecision!="Rejected"])
fNamesAd<-gsub("`","",fNamesAd)
g1<-intersect(fNamesAd,Genes)

library(VennDiagram)   
library(gplots)
col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6",
        "#E2D200","#B40F20","#273046","#FD6467")
col2=c("#5BBCD6","#B40F20")
col3=c("#5BBCD6","#B40F20","#E2D200")
col5=c("#5BBCD6","#FD6467","#E2D200","#00A08A","#F98400")
barplot(1:5,col=col5)

#'Trajectory Inference' =Genes,
#"Boruta selection"=fNamesAd
venn.diagram(
  x = list('scVelo' = velo$Gene,
           'Monocle3' = Monocle3$gene_short_name,
           'Monocle2' = Monocle2$x,
           'slingshot' =slingshot$x,
           'destiny' = destiny$x
  ),
  filename = 'VN-features.jpg',
  col = "black",
  fill = col5,
  alpha = 0.5,
  cex = 1.2,
  cat.col = col5,
  cat.cex = 2,
  #main = "gene selection",
  #main.cex = 1.2,
  label.col = "white",
  fontfamily = "serif",
  # fontface = "bold",
  cat.fontfamily = "serif",
  # cat.fontface = "bold",
  margin = 0.18
  # cat.dist = c(0.03, 0.03),
  # cat.pos = c(-20, 20)
)

venn.diagram(
  x = list('Monocle3' = Monocle3$gene_short_name,
           'slingshot' =slingshot$x),
  filename = 'VN-TI.png',
  col = "black",
  fill = col2,
  alpha = 0.5,
  cex = 1.5,
  cat.col = 'black',
  cat.cex = 0,
  cat.fontface = "bold",
  margin = 0.1,
  main = "gene selection",
  main.cex = 1.2
)
union()
library(eulerr)
col3=c("#fbb4ae", "#b3cde3", "#ccebc5")

vd <- euler(c(A= 116,
              B=1181,
              C=0, 
              "A&B" = 80, "A&C" = 7, "B&C" = 362,
              "A&B&C" = 77))
vd <- euler(c(A= 123,
              B=601,
              C=717, 
              "A&B" = 99, "A&C" = 23, "B&C" = 225,
              "A&B&C" = 35))
jpeg(file="VN-features-2.jpeg",width=6,height=6,units = "in", res = 1000)
plot(vd,
     fills = list(fill = col3, alpha = 0.6),
     labels = list(col = "white", font = 4), 
     edges = T,
     quantities = TRUE)
dev.off()

vd <- euler(c(A= 740,
              B=700,
              "A&B" = 260))

jpeg(file="VN-TI.jpeg",width=4,height=4,units = "in", res = 500)
plot(vd,
     fills = list(fill = col2, alpha = 0.5),
     labels = list(col = "white", font = 4), 
     edges = FALSE,
     quantities = TRUE)
dev.off()

library(UpSetR)
all<-union(Genes)
df<-data.frame(row.names = Genes,
               features=ifelse(Genes%in%fNamesAd,1,0),
               slingshot=ifelse(Genes%in%slingshot$x,1,0),
               Monocle3=ifelse(Genes%in%Monocle3$gene_short_name,1,0),
               scVelo=ifelse(Genes%in%velo$Gene,1,0),
               Monocle2 = ifelse(Genes%in%Monocle2$x,1,0),
               destiny = ifelse(Genes%in%destiny$x,1,0))
which(is.na(Genes))
which(is.na(velo$Gene)) 
?upset()
colnames(df)
col8 <- c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
          "Myeloid-CCL5"="#B40F20","Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
          "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A")
setcol<-c("features"="#FF0000","slingshot"="#E2D200","Monocle3"="#5BBCD6",'destiny' = "#FD6467","scVelo"="#00A08A",'Monocle2' = "#F98400")
jpeg(file="VN-features-2.jpeg",width=8,height=4,units = "in", res = 1000)
eps(file="VN-features-2.eps",width=8,height=4) 
setEPS()
postscript("VN-features-2.eps")
upset(df, nsets = 6, sets = c("features","Monocle3",'Monocle2',"slingshot",'destiny',"scVelo"),
      mb.ratio = c(0.6,0.4),#调整上下两部分的比例
      order.by = c("degree"),keep.order=T,
      #nintersects=10,
      sets.bar.color = setcol,point.size=4,line.size=1,text.scale = 1.5)
# queries = list(list(query = intersects,
#           params = list("features","Monocle3","slingshot","scVelo"), 
#           color ="#FF0000", active = T),
#           list(query = intersects,
#                params = list("features","Monocle3","slingshot"), 
#                color ="#FF0000", active = T),
#           list(query = intersects,
#                params = list("features","Monocle3"),
#                color ="#FF0000", active = T),
#           list(query = intersects,
#                params = list("features","Monocle3","scVelo"),
#                color ="#FF0000", active = T),
#           list(query = intersects,
#                params = list("features","slingshot"),
#                color ="#FF0000", active = T),
#           list(query = intersects,
#                params = list("features","scVelo"), 
#                color ="#FF0000", active = T)))#为按频率排序
dev.off()
#################cor########
jpeg(file="corTI-allcell.2.10.jpeg",width=7,height=5.5,units = "in", res = 1000)
corrplot(cordfa, type = "upper", tl.pos = "tp",
         tl.cex=1.5,col = colorRampPalette(c("#46ACC8","white","#B40F20"))(50),
         cl.cex = 0.8,#order = "hclust", addrect = 2,
         tl.col = "black",family="serif")
corrplot(cordfa, add = TRUE, type = "lower", method = "number",
         col = colorRampPalette(c("#46ACC8","white","#B40F20"))(50), 
         cl.cex = 4,family="serif",
         diag = FALSE, tl.pos = "n", cl.pos = "n")
dev.off()

load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
Mono_Macr<-subset(x = Mono_Macr,subset = ann!="Macr-CCL5")
matrix<-Mono_Macr@assays$RNA@data
matrix[1:5,1:5]
spp1<-matrix["SPP1",]
mif<-matrix["MIF",]
cd74<-matrix["CD74",]

cor(spp1,mif)#0.1999671
cor(spp1,cd74)
cor(cd74,mif)
scdata=data.frame(MIF=mif,SPP1=spp1,CD74=cd74,row.names = colnames(matrix))
rm(matrix,Mono_Macr)

exprMatrix = read.table(file = "D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\prognosis\\seq_all_symbol_tpm.txt",
                        header=TRUE,row.names=1, as.is=TRUE)
exprMatrix = read.table(file = "D:\\Bioinfrolf\\Bioinfrolf\\valid\\HCCDB\\merge\\merge_symbol_tpm.txt",header=TRUE,row.names=1, as.is=TRUE)
Genes<-read.table("D:/Bioinfrolf/Bioinfrolf/GSEA/impMatrix_201.txt",header = T)
index=intersect(rownames(exprMatrix),Genes$Feature)

ex=exprMatrix[index,]
ex["SPP1",]->spp1
ex["MIF",]->mif
ex["CD74",]->cd74
seqdata<-rbind(SPP1=spp1,MIF=mif,CD74=cd74)%>%t()%>%as.data.frame()
arrydata<-rbind(SPP1=spp1,MIF=mif,CD74=cd74)%>%t()%>%as.data.frame()
rm(exprMatrix,ex)

cor(ndata[,1:2])
cor(as.numeric(spp1),as.numeric(mif))#0.3233247
ce<-cor(t(ex))

data<-scdata
x=19
OPT<-function(x,data){
  ndata<-c(0,0,0)%>%as.data.frame()
  data<-data[order(data[,1],decreasing = T),]
  data<-data[data$MIF!=0,]
  #data<-data[!c(data$MIF==0 &data$SPP1==0),]
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

exprMatrix.2<-exprMatrix.2[,51:ncol(exprMatrix.2)]
data<-exprMatrix.2[c("CXCL10","CXCL9","CXCL10"),]%>%t()%>%as.data.frame()
rt<-as.data.frame(results0[51:nrow(results0),])
rownames(rt)<-make.names(rownames(rt))
cdt<-cor(data$CXCL10,rt$`T cells CD8`)
rt1<-rt[rt$`P-value`<0.05,]
seqdata3<-rbind(CXCL10=immuneScore$CXCL10,CD8=score2,NK=score3)%>%t()%>%as.data.frame()
seqdata3<-seqdata3[order(seqdata3$CXCL10,decreasing = T),]
seqdata10<-rbind(CXCL10=immuneScore$CXCL10,CD8=score2,NK=score3)%>%t()%>%as.data.frame()
seqdata10<-seqdata10[order(seqdata10$CXCL10,decreasing = T),]
seqdata9<-rbind(CXCL9=immuneScore$CXCL9,CD8=score2,NK=score3)%>%t()%>%as.data.frame()
seqdata9<-seqdata9[order(seqdata9$CXCL9,decreasing = T),]

seqcdt3<-OPTy(30,seqdata3)
seqcdt10<-OPTy(30,seqdata10)
seqcdt9<-OPTy(30,seqdata9)
arrycdt<-OPTy(30,arrydata)
sccdt<-OPTy(20,scdata)
scdata.op<-OPT(19,scdata)

save(seqcdt,arrycdt,sccdt,file = "cor.data")

jpeg(file="CXCL9-NK.jpeg",width=8,height=6,units = "in", res = 2000)
ggplot(data=seqdata3, aes(x=CXCL9, y=score3)) +
  geom_point(alpha=0.4, size=2) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("CXCL9") + ylab("CD8 T cell") + 
  scale_colour_manual(values = c("darkgreen", "grey", "brown4"),limits = c("DOWN","NOT","UP"))+
  theme_bw()+theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_smooth(aes(x=CXCL9, y=score3),method = lm,linetype=1,se=T,span=5)+
  annotate("text",x=10,y=0.15,label="R:0.746",parse=T)+xlim(c(0,250))+ylim(c(0,0.12))
dev.off()

ggplot(data = data.demon2,aes(x=endpoint,y=event.rate))+##数据列
  geom_point(aes(size=event.rate,color=IN.C),alpha=0.7,show.legend = TRUE)+##点图大小、颜色、透明度、图例
  geom_smooth(aes(x=endpoint,y=event.rate,color=IN.C),method = lm,linetype=1,se=FALSE,span=1)+##趋势线、颜色、方法、线型、置信区间
  guides(color=guide_legend(title=NULL))+##去除图例标题
  scale_size(range = c(1, 10),guide=FALSE)+##气泡大小区间，图例标题去除
  scale_color_manual(values=c("CornflowerBlue","Gold"),
                     breaks=c("INevent.rate","Cevent.rate"),
                     labels=c("1","2"))+##气泡颜色、图例名称修改)+##气泡颜色
  labs(x='Endpoint',y='Increase ratio')+##横纵坐标
  geom_text_repel(aes(label=event.rate))+##气泡标注
  annotate("text",x=1.0,y=2.8,label="atop(Y==3*X+5,R^2==0.9)",parse=T)+##添加公式
  #theme(legend.title = element_blank())+
  theme_bw();p.demon##设置背景theme

jpeg(file="Seq-genecor.jpeg",width=8,height=6,units = "in", res = 1000)
hseq<-pheatmap::pheatmap(ce,#scale = "row",
                         #number_color="gray",
                         # display_numbers = matrix(ifelse(ce < 0.001, "***", ifelse(
                         #   matpval < 0.01,"**",ifelse(
                         #     matpval < 0.05,"*","-"
                         #   )
                         # )),nrow(matBP)),
                         show_colnames =F,
                         show_rownames = T,
                         cluster_row = T, cluster_col = T,
                         border_color ='white',
                         angle_col=45,fontsize=10,fontsize_col=10,fontsize_row=6,
                         color = colorRampPalette(rev(brewer.pal(n = 11, name =
                                                                   "RdYlBu")))(100),
                         main ="Seq genecor") 
dev.off()
h$tree_row$labels
h$kmeans
ce["SPP1","MIF"]
hseq<-h
save(hseq,hary,file = "cor.Rdata")
#################accuracy##############################
library(tidyr)
library(tibble)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)
library(gg.gap)
library(ggprism)
library(ggpubr)

rm(list = ls())
setwd("D:/Bioinfrolf/HCC-SC/Myeloid/model")
databar=read.csv(file='accr-1.csv',header = T,stringsAsFactors = F)
databar$pair<-rep(1:3,6)
col2<-c("#5BBCD6","#B40F20")#c("#FF0000","#00A08A")

q1<-ggplot(data=databar, mapping=aes(x = sample, y = Accuracy,
                                     fill=sample))+
  geom_bar(stat="identity",position=position_dodge(0.75),
           show.legend = T, width = 0.6)+
  scale_fill_manual(values = col2)+
  geom_point(size=3,aes(x = sample, y = Accuracy))+
  geom_line(aes(group = pair), lwd = 0.5)+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 10, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),)+
  scale_y_continuous(labels=percent,breaks =seq(0,1,0.04) )
q1

q2=q1+facet_wrap(~method,as.table=F)+theme(
  strip.background = element_rect(fill = 'white',color = 'white'),
  #strip.placement = 'outside',
  strip.text = element_text(color = 'black',size = 15)
)
q2

p2 = gg.gap(plot = q2,
            segments = c(0.1, 0.8),
            tick_width = 0.05,
            rel_heights = c(0.2, 0, 0.3),# 设置分隔为的三个部分的宽度
            ylim = c(0, 1)
)


jpeg(file="accuracy-1.jpeg",width =6,height = 4,units = "in", res = 2000)
p2
dev.off()

#################violin################
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
Mono_Macr<-subset(x = Mono_Macr,subset = ann!="Macr-CCL5")
setwd("D:/Bioinfrolf/Bioinfrolf/GSEA/")
immunomodulator<-read.table("D:/Bioinfrolf/Bioinfrolf/GSEA/immunomodulator.txt",header = F)
colnames(immunomodulator)<-c("protein","type","gene")

sce<-Mono_Macr
vio <-as.matrix(sce@assays$RNA@data)
vio<-t(vio)
vio<-as.data.frame(vio)
immunomodulator<-immunomodulator[immunomodulator$gene%in%colnames(vio),]
chemokine<-immunomodulator$gene[immunomodulator$type=="chemokine"]
receptor<-immunomodulator$gene[immunomodulator$type=="receptor"]
MHC<-immunomodulator$gene[immunomodulator$type=="MHC"]
Immunoinhibitor<-immunomodulator$gene[immunomodulator$type=="Immunoinhibitor"]
Immunostimulator<-immunomodulator$gene[immunomodulator$type=="Immunostimulator"]

features<-intersect(colnames(vio),M1)#intersect(colnames(vio),c(res,inf))
features<-c("CD44","CXCR4","CD74","MIF")
vio<-vio[,features]
# s<-colSums(vio)
# which(s=)
# apply(vio, 2, max)-apply(vio, 2, min)
#vio<-as.matrix(vio)
# colnames(vio)<-immunomodulator$type
# rownames(vio)<-Mono_Macr$ann
vio<-cbind(sce$ann,vio)
library(reshape2)
new_vio <- melt(vio)
colnames(new_vio)<-c("celltype","marker","value")
head(new_vio)

col8 <- c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
          "DC-CCL5"="#B40F20","Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
          "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A")
annocol2=c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
           "Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
           "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A",
           "MDM-APOE/APOC1" ="#B40F20","MDM-CD4"="#E6A0C4","MDM-VCAN"="#354823",
           "MDM-HLA-DQA1"="#F2AD00","MDM-THBS1"="#C6CDF7","KC-VCAM1"="#238443")

jpeg(file="vio-kc.jpeg",width=6,height=4,units = "in", res = 1000)
p<-ggplot(new_vio,aes(x = celltype,y = value)) +
  geom_violin(aes(fill = celltype,color = celltype),show.legend = F) +
  geom_boxplot(width=0.07, color="gray", alpha=0)+
  theme_bw() +
  coord_flip() +
  facet_grid(~marker,scales = 'free') +
  scale_x_discrete(position = 'top') +
  xlab('') + ylab('') +
  scale_fill_manual(values=annocol2) +scale_color_manual(values = rep(annocol2))+
  theme(panel.grid = element_blank(),
        # 分面x标签背景
        strip.background.x = element_blank(),
        panel.border = element_rect(size = 0.5),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(angle = 90,vjust = 0.25,hjust = 0,size = 8),
        # 分面panel x轴上间距
        panel.spacing.x = unit(0,'mm')
  )
dev.off()
library(grid)
jpeg(file="vio-MIF.jpeg",width=4,height=4,units = "in", res = 1000)
p#print(p, vp=viewport(angle=90))
dev.off()

vio<-data.frame(value=sce@assays$RNA@data["MIF",],celltype=sce$ann,
                row.names = colnames(sce@assays$RNA@data))
group_by(vio,celltype)
df = tapply(vio$value, INDEX=as.factor(as.character(vio$celltype)), FUN=mean)%>%as.data.frame()
colnames(df)<-"value"
df$celltype<-rownames(df)
df<-df[c(7,6,3,1,4,5,2),]
df$celltype<-factor(df$celltype,levels = df$celltype)
vio$celltype<-factor(vio$celltype,levels = df$celltype)

jpeg(file="bar-MIF.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(df,aes(x = celltype,y = value,fill=celltype)) +
  geom_bar(stat="identity",show.legend = F,width = 0.6) +
  xlab('') + ylab('MIF Average Expression') +
  scale_fill_manual(values=annocol2) +scale_color_manual(values = rep(annocol2))+
  #geom_point(size = 2) +  #绘制样本点
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 10),
        axis.text = element_text(size = 10, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size = 12))#+ylim(0,0.1)
dev.off()

jpeg(file="vio-MIF1.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(vio,aes(x = celltype,y = value)) +
  geom_violin(aes(fill = celltype,color = celltype),show.legend = F) +
  xlab('') + ylab('MIF Expression') +
  scale_fill_manual(values=annocol2) +scale_color_manual(values = rep(annocol2))+
  geom_boxplot(width=0.1, color="black", alpha=0)+
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 10),
        axis.text = element_text(size = 10, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1,size = 12))
dev.off()
