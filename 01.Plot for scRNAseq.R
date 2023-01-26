rm(list=ls())
library(ggplot2)
library(ggprism)
library(ggsci)
library(patchwork)
library(reshape2)
library(Seurat)
library(tidyverse)
library(patchwork)
library(pheatmap)
library("RColorBrewer")
library(wesanderson)
options(stringsAsFactors = F)
setwd("D:/Bioinfrolf/new figs/")
#########################color###############################
brewer.pal.info
display.brewer.all()
cols<-brewer.pal(11,"PRGn")
cols<-cols[c(2,10)]
brewer.pal(9,"Set1")
annocol2=c("Macr-APOE" ="#FF0000","Macr-RPs"="#F98400","Mono-VCAN"="#273046","Mono-CD16"="#FD6467",
           "Macr-STAT1"="#E2D200","Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A",
           "hMDM-APOE" ="#B40F20","hMDM-CD14"="#E6A0C4","hMDM-VCAN"="#354823",
           "hMDM-HLA-DQA1"="#F2AD00","hMDM-THBS1"="#C6CDF7","hKC-VCAM1"="#238443")
annocol1<-c("MDM-APOE" ="#B40F20","MDM-CD14"="#E6A0C4","MDM-VCAN"="#354823",
            "MDM-HLA-DQA1"="#F2AD00","MDM-THBS1"="#C6CDF7","KC-VCAM1"="#238443")
ann_color<-list(Tissue=c("HN1"="green", "Normal"="#00A08A","TumorEdge"="#F98400" ,"TumorCore"="#FF0000"),
                CellType=annocol2)
col2<-c("#B40F20","#5BBCD6")
col3<-c("AN"="#00A08A","PT"="#F98400","CT"="#FF0000")

col4<-c("Hepatocytes"="#E6A0C4","Fibroblasts"="#F2AD00","Myeloid"="#B40F20","Endothelial Cells"="#7294D4")
col6<-c(wes_palettes$FantasticFox1[c(2,5)],wes_palettes$Darjeeling1[c(1,2,4,5)])
col8 <- c("Macr-APOE" ="#FF0000","Macr-RPs"="#F98400","Mono-VCAN"="#273046",
         "Mono-CD16"="#FD6467","Macr-STAT1"="#E2D200",
          "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A")
col12<-c("NK cells"="#5BBCD6","CD8 T cells"="#00A08A",
         "Myeloid cells"="#FF0000","ILCs"="#B40F20","Mast cells"="#7294D4",
         "Tregs"="#F98400","CD4 T cells"="#F2AD00","HSPs CD4 T cells"="#E2D200",
         "HSPs T/NK cells"="#C6CDF7","T-NK-Cycle"="#273046",
         "B cells"="#FD6467","Plasma cells"="#E6A0C4")
col14.1<-c("NK-GNLY"="#5BBCD6","T-CD8-GZMK"="#00A08A","Myeloid-S100A9"="#FF0000","T-CD4-IL7R"="#F2AD00",
         "T/NK-Cycle-STMN1"="#273046","B-CD79A"="#FD6467","Myeloid-APOC1"="#B40F20","Mast-TPSB2"="#7294D4",
         "B-Plasma-MZB1"="#E6A0C4","Treg-CTLA4"="#F98400","T-CD4-HSPA6"="#E2D200","NK-TM4SF1"="#D8A499",
         "ILC-LILRA4"="#354823","T/NK-HSPA6"="#C6CDF7")
col13 <- c("Macr-APOE" ="#FF0000","Macr-TMSB4X"="#F98400","Mono-VCAN"="#273046","Mono-CD16"="#FD6467",
           "Macr-STAT1"="#E2D200","Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A","Myeloid-CCL5"="#B40F20",
           "DC-CD1C"="#D8A499","DC-CLEC9A"="#E6A0C4", "DC-MKI67"="#C6CDF7",
           "DC-LAMP3"="#7294D4", "pDC-IGKC"="#354823")
col14 <- c("Macr-APOE" ="#FF0000","Macr-RPs"="#F98400","Mono-VCAN"="#273046","Mono-CD16"="#FD6467",
          "Macr-STAT1"="#E2D200","Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A","Myeloid-CCL5"="#B40F20",
          "DC-CD1C"="#D8A499", "DC-CD1C-RPs"="#F2AD00", "DC-CLEC9A"="#E6A0C4", "DC-MKI67"="#C6CDF7",
          "DC-LAMP3"="#7294D4", "pDC-IGKC"="#354823")

col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6",
        "#E2D200","#B40F20","#273046","#FD6467")
col14<-c(col8,c("#354823","#E6A0C4","#F2AD00","#C6CDF7","#D8A499","#7294D4"))

#########################scRNAseq QC############################
jpeg(file="umap_GSE-hn.jpg",width=10,height=4,units = "in", res = 1000)
DimPlot(sce, reduction = "tsne", group.by = "GSE", pt.size =0.5, split.by = 'GSE',cols = c("#B40F20","#E2D200","#5BBCD6"),)
dev.off()

sce$ann<-factor(sce$ann,levels = c("Mono-VCAN","Mono-CD16","Macr-SPP1","Macr-STAT1","Macr-VCAM1","Macr-APOE",
                                   "Macr-RPs","DC-CD1C-RPs","DC-CD1C","DC-CLEC9A","DC-MKI67","DC-LAMP3",
                                   "pDC-IGKC","Myeloid-CCL5"))
meta<-sce@meta.data[,c("nFeature_RNA", "nCount_RNA","percent.ribo","ann")]#"percent.mt","percent.ribo"
head(meta)
jpeg(file="myeloid-nFeature.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(meta, aes(x = ann, y = nFeature_RNA)) + 
  geom_violin(aes(fill = ann,color = ann),show.legend = T) +
  geom_boxplot(width=0.5, color="black", alpha=0)+
  labs(title ="HCC Myeloid Cells",x = '',y = "nFeature")+
  theme_bw() +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
  )+scale_fill_manual(values=col14)+scale_color_manual(values=col14)+ylim(1,5000)
dev.off()
cols<-annocol1
jpeg(file="hn-nCount.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(meta, aes(x = ann1, y = nCount_RNA)) + 
  geom_violin(aes(fill = ann1,color = ann1),show.legend = T) +
  geom_boxplot(width=0.5, color="black",alpha=0)+
  labs(title ="HCC Myeloid Cells",x = '',y = "nCount")+
  theme_bw() +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
  )+scale_fill_manual(values=cols)+scale_color_manual(values=cols)#+ylim(1,20000)
dev.off()

jpeg(file="vio-ribo-1.jpeg",width=4,height=4,units = "in", res = 1000)
ggplot(meta, aes(x = ann, y = percent.ribo/100)) + 
  geom_violin(aes(fill = ann,color = ann),show.legend = T) +
  geom_boxplot(width=0.5, color="black", alpha=0)+
  labs(title ="HCC MNPs",x = '',y = "percent.ribo")+
  theme_bw() +coord_flip() +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
  )+scale_fill_manual(values=col14)+scale_color_manual(values=col14)+#+ylim(1,5000)
  scale_y_continuous(labels = scales::label_percent())+NoLegend()
dev.off()

#########################scRNAseq celltype bar plot####################
head(sce@meta.data)
table(sce$GSE)
metadata<-data.frame(cell=rownames(sce@meta.data),celltype=sce$ann,GSE=sce$GSE)
rownames(metadata)<-metadata$cell

dat<-table(metadata$GSE,metadata$celltype)%>%as.data.frame()
colnames(dat)<-c("GSE","celltype","percent")
dat$percent<-ifelse(dat$GSE=="GSE115469",dat$percent/1127,ifelse(dat$GSE=="GSE136103",dat$percent/1841,dat$percent/71))
dat$percent<-ifelse(dat$GSE=="GSE156625",dat$percent/5470,dat$percent/4911)
dat$celltype<-factor(dat$celltype,levels = c("Mono-VCAN","Mono-FCGR3A","Macr-VCAM1","Macr-CXCL10",
                                             "Macr-APOE","Macr-APOC1","Macr-SPP1"))

current.cluster.ids <-levels(as.factor(sce$tissue_sub))
new.cluster.ids <-c("AN","CT","PT")
sce$tissue_sub <- plyr::mapvalues(x = sce$tissue_sub,
                                           from = current.cluster.ids, to = new.cluster.ids)
metadata<-data.frame(cell=rownames(sce@meta.data),celltype=sce$global,tissue=sce$tissue_sub)
rownames(metadata)<-metadata$cell
dat<-table(metadata$tissue,metadata$celltype)%>%as.data.frame()
colnames(dat)<-c("tissue","celltype","percent")
dat$tissue<-factor(dat$tissue,levels = c("AN","PT","CT"))
dat$celltype<-factor(dat$celltype,levels = c("Tregs","HSPs T/NK cells","HSPs CD4 T cells","T-NK-Cycle",
                                             "CD4 T cells","CD8 T cells","NK cells",
                                         "B cells","Plasma cells","Myeloid cells","ILCs","Mast cells"))
jpeg(file="bar-myeloid-GSE.jpeg",width=6,height=4,units = "in", res = 1000)
ggplot(dat,aes(x = celltype,y = percent)) +
  geom_bar(aes(fill = GSE),
           stat = 'identity',position = 'fill',
           show.legend = T) +
  theme_bw() +
  theme(axis.line = element_line(colour = 'black', size = 0.5), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.y= element_text(size = 12),
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
  )+scale_y_continuous(labels = scales::label_percent())+
  labs(title ="Liver MNPs",x = '',y = "percent")+scale_fill_manual(values=col2)
dev.off()


#########################scRNAseq UMAPplot####################
rm(list = ls())
table(sce$ann1)
if(cancer){
  umap<-rbind(normal@reductions$umap@cell.embeddings,
              para@reductions$umap@cell.embeddings,
              core@reductions$umap@cell.embeddings)%>%
    as.data.frame()%>% cbind(tx = c(normal$PNC,para$PNC,core$PNC))
  new.cluster.ids <- c("CT","AN","PT") 
  current.cluster.ids<-levels(as.factor(umap$tx))
  umap$tx<- plyr::mapvalues(x =umap$tx,from = current.cluster.ids, to = new.cluster.ids)
  table(umap$tx)
}
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/04.HN_HCC.Rdata")
{
  table(sce$ann1)
  new.cluster.ids <- c("hMDM-APOE","hMDM-CD14","hMDM-HLA-DQA1","hMDM-THBS1","hMDM-VCAN","hKC-VCAM1",
                       "Macr-APOE","Macr-RPs","Macr-SPP1","Macr-STAT1","Macr-VCAM1", 
                       "Mono-CD16","Mono-VCAN") 
  
  current.cluster.ids<-levels(as.factor(sce$ann1))
  sce$ann1<- plyr::mapvalues(x =sce$ann1,from = current.cluster.ids, to = new.cluster.ids)
  save(sce,file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/04.HN_HCC.Rdata")
  sce$sample%>%table()
}
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/01.ALLimmune.combined.Rdata")
if(F){
  rm(immune.combined)
  head(sce@meta.data)
  table(sce$ann)
  table(sce$global)
  hn$ann1%>%unique()
  sce$ann1%>%unique()
  umap = sce@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = sce$cluster)
  cols<-col4#annocol2#col12#
  p<-ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
    geom_point(size = 0.2, alpha = 1,show.legend = F)+
    labs(title ="cells for CellChat",x = 'UMAP1',y = "UMAP2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)+NoLegend()
    )+scale_color_manual(values=col3)
  
  jpeg(file = "umap-cancer-tissue.jpg", width = 4.5, height = 4, units = "in", res =1000)
  p#LabelClusters(plot = p, id = "tx",size = 4,color = "black",repel = T)
  dev.off()

  umap = sce@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = sce@meta.data$ann)
  jpeg(file = "umap-MNPs.jpg", width = 2.75, height = 3, units = "in", res =1000)
  ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
    geom_point(size = 0.2, alpha = 1,show.legend = F)+
    labs(title ="cells for CellChat",x = 'UMAP1',y = "UMAP2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 15, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)+NoLegend()
    )+scale_color_manual(values=annocol2)
  dev.off()
  
  jpeg(file = "UMAP-MNPs-1.jpg", width = 6, height = 4, units = "in", res =1000)
  ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
    geom_point(size = 0.2, alpha = 1,show.legend = T)+
    labs(title ="HCC MNPs",x = 'UMAP1',y = "UMAP2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)
    )+scale_color_manual(values=annocol2)
  dev.off()
  
  tsne = sce@reductions$tsne@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = sce@meta.data$ann1)
  jpeg(file = "HN-tsne-MNPs.jpg", width = 6, height = 4, units = "in", res =1000)
  ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = tx)) + 
    geom_point(size = 1, alpha = 1,show.legend = T)+
    labs(title ="Healthy Liver cells",x = 'tSNE1',y = "tSNE2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)
    )+scale_color_manual(values=annocol1)
  dev.off()
  
  umap = sce@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = sce@meta.data$SingleRCluster)
  jpeg(file = "umap-myeloids-singleR.jpg", width = 6, height = 4, units = "in", res = 1000)
  ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
    geom_point(size = 0.2, alpha = 1)+#,key_glyph="polygon"
    labs(title ="The Blueprint/ENCODE reference", x = 'UMAP1',y = "UMAP2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)
    )+scale_color_manual(values=col3)
  dev.off()
  
  umap = immune.combined@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = immune.combined$GSE)
  umap1 = umap[umap$tx=="GSE140228",]
  umap2 = umap[umap$tx!="GSE140228",]
  jpeg(file = "umap-immune-GSE-s.jpg", width = 6, height = 4, units = "in", res = 1000)
  ggplot(umap2, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
    geom_point(size = 0.01, alpha = 1)+#,key_glyph="polygon"
    labs(title ="GSE156625", x = 'UMAP1',y = "UMAP2")+
    theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.title = element_blank(), 
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12)
    )+scale_color_manual(values="#5BBCD6")
  dev.off()
}

#########################scRNAseq Feature plot#####################
features=list("CD8 T cells"=c("CD8A","GZMK"),"NK cells"=c("GNLY","KLRD1"),
              "HSPs T/NK cells"=c("IFNG","HSPA6"),"HSPs CD4 T cells"= c("DNAJB1","CD69"),
              #"exhausted T/NK cells"=c("HAVCR2","TIGIT"),"exhausted CD4 T cells"= c('NR4A1',"CD69"),
              "Tregs"=c("CTLA4","FOXP3"),"T-NK-Cycle"=c("MKI67","STMN1"),"Plasma cells"=c("MZB1","CD38"),
              "B cells"=c("CD79A","CD24"),"Mast cells"=c("TPSB2","KIT"),
              "Myeloid cells"=c("LYZ","CST3"),"ILCs"=c("LILRA4","IRF4"))
levels(sce)
levels(sce$global)
sce$global<-factor(sce$global,levels = c("CD8 T cells","NK cells","HSPs T/NK cells","HSPs CD4 T cells","CD4 T cells",
                                         "Tregs","T-NK-Cycle","Plasma cells","B cells","Mast cells","Myeloid cells","ILCs"))

jpeg(file="Dotplot-immune.jpeg",width =9.5,height = 4,units = "in", res = 1000)
DOTplot(sce, features = features,group.by="global",assay ="integrated",
        cols = c("#5BBCD6","#B40F20")) 
dev.off()

features=list("Mono-VCAN"=c("VCAN","S100A9"),"Mono-CD16"=c("FCGR3A","LILRA5"),"Macr-SPP1"=c("SPP1","NUPR1"),
              "Macr-STAT1"=c("STAT1","CXCL10"),"Macr-VCAM1"=c("VCAM1","CD5L"),"Macr-APOE"=c("APOE","TREM2"),
              "Macr-RPs"=c("ATP5E","RPLP1"),"DC-CD1C-RPs"=c("RPS8","HLA-DPB1"),"DC-CD1C"=c("CD1C","CD1E"),
              "DC-CLEC9A"=c("CLEC9A","DNASE1L3"),"DC-MKI67"=c("MKI67","STMN1"),"DC-LAMP3"=c("LAMP3","CCR7"),
              "pDC-IGKC"=c("IGKC","GZMB"),"Myeloid-CCL5"=c("CCL5","CD3D"))
sce$ann<-factor(sce$ann,levels = c("Mono-VCAN","Mono-CD16","Macr-SPP1",
                                   "Macr-STAT1","Macr-VCAM1","Macr-APOE",
                                   "Macr-RPs","DC-CD1C-RPs","DC-CD1C",
                                   "DC-CLEC9A","DC-MKI67","DC-LAMP3","pDC-IGKC","Myeloid-CCL5"))
features<-list("Kupffer cells"=c("VCAM1","MARCO","VSIG4"),
            "Non-inflammatory cells"=c("CD163","MRC1","FCGR3A"),
            "Inflammatory cells"=c("S100A8","S100A9","S100A6","S100A12","VCAN" ,"LYZ","FCN1"),
            "Monocyte-derived cells"=c("CCR2","ITGAM","CD14"))
features<-list("a"=c("VCAM1","MARCO","VSIG4"),
               "b"=c("CD163","MRC1","FCGR3A"),
               "c"=c("S100A8","S100A9","S100A6","S100A12","VCAN" ,"LYZ","FCN1"),
               "d"=c("CCR2","ITGAM","CD14"))

jpeg(file="Dotplot-HN.jpeg",width =8,height = 2.75,units = "in", res = 1000)
DOTplot(sce, features = features,group.by="ann1",assay ="RNA",
        cols = c("#5BBCD6","#B40F20"))
dev.off()

jpeg(file="FEATUREplot-TAMs-1.jpeg",width =10,height = 8,units = "in", res = 1000)
FEATUREplot(sce, features = c(
  "C1QA","C1QB","C1QC","MRC1","CD163","MS4A4A","MIF","TREM2","SPP1"),
  cols = c("light grey","#B40F20"), min.cutoff = "q3")
dev.off()

features<-c("CXCL2","CXCL3", "CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL16")
jpeg(file="Dotplot-CXCLs.jpeg",width =5.5,height = 2.75,units = "in", res = 1000)
DOTplot(sce, features = features,group.by="ann",assay ="RNA",
        cols = c("#5BBCD6","#B40F20"))
dev.off()

#########################scRNAseq Feature violin###############
features=c("HAVCR2","TIGIT",'NR4A1',"CD160","PDCD1","EOMES","TBR2","LAG3","FDC")
library(reshape2)
unique(immune.combined$global)
cells<-rownames(immune.combined@meta.data[immune.combined$global%in%c("CD8 T cells","NK cells","HSPs T/NK cells",
                                                                      "HSPs CD4 T cells","CD4 T cells",
                                                                      "Tregs", "T-NK-Cycle"),])
vio <-immune.combined@assays$RNA@data[,cells]
vio<-vio[intersect(features,rownames(vio)),]
vio<-t(as.matrix(vio))
vio<-as.data.frame(vio)
vio<-cbind(vio,immune.combined@meta.data[cells,"global"])
new_vio <- melt(vio)
colnames(new_vio)<-c("celltype","marker","value")
head(new_vio)

jpeg(file="vio.jpeg",width=12,height=8,units = "in", res = 2000)
ggplot(new_vio,aes(x = celltype,y = value)) +
  geom_violin(aes(fill = celltype,color = celltype),show.legend = F) +
  geom_boxplot(width=0.5, color="black", alpha=0)+
  theme_bw() +
  coord_flip() +
  facet_grid(~marker,scales = 'free') +
  scale_x_discrete(position = 'top') +
  xlab('') + ylab('') +
  scale_fill_manual(values=col12) +scale_color_manual(values = rep(col12))+
  theme(panel.grid = element_blank(),
        # 分面x标签背景
        strip.background.x = element_blank(),
        panel.border = element_rect(size = 0.5),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = 'bold',size = 5),
        # 分面x标签旋转角度
        strip.text.x = element_text(angle = 45,vjust = 0.25,hjust = 0,face = 'bold',size = 4),
        # 分面panel x轴上间距
        panel.spacing.x = unit(0,'mm')
  )
dev.off()
#########################scRNAseq tissue pie################
table(sce$sample)
hn<-subset(sce,sample=="HN1")

info = table(hn$ann1)%>%as.numeric()
names = hn$ann1%>%as.factor()%>%levels()
names(info)<-names
cols = annocol2
piepercent = paste(round(100*info/sum(info)), "%")

jpeg(file="HN-pie-tissue.jpeg",width =10,height =5,units = "in", res = 1000)
pie(info, labels=piepercent, main = "HN", col=cols,cex=2, border = NA)
legend("topright", names, cex=1.5, fill=annocol2,border =NA,bty="n")
dev.off()

info = c(71,1127,1841)
names = c("GSE129933","GSE115469","GSE136103")
cols = c("#E2D200","#B40F20","#5BBCD6")
piepercent = paste(round(100*info/sum(info)), "%")
jpeg(file="hn-pie-GSE.jpeg",width =8,height =5,units = "in", res = 1000)
pie(info, labels=piepercent, main = "data source", col=cols,cex=2, border = NA)
legend("topright", names, cex=1.5, fill=cols,border =NA,bty="n")
dev.off()


head(sce@meta.data)
table(hn$ann)
table(hn$GSE)
head(sce@meta.data)
sce[,sce$GSE=="GSE156625"]
sce[,sce$GSE=="GSE140228"]


save(list = ls(),file = "hnsce.Rdata")

metahn<-sce@meta.data[,5:10]
table(sce$ann)
table(sce$GSE)
table(sce$tissue_sub)
table(sce$patient)
table(sce$tissue)
info = c(3654,3262,1209,62)
names = c("GSE156625","GSE140228","GSE115469","GSE129933")
cols = c("#FF0000","#00A08A", "#FAD510","#F98400")#c("#FF0000","#00A08A")
piepercent = paste(round(100*info/sum(info)), "%")

jpeg(file="pie-GSE-hn.jpeg",width =8,height =5,units = "in", res = 2000)
pie(info, labels=piepercent, main = "data source", col=cols,cex=2, border = NA)
legend("topright", names, cex=1.5, fill=cols,border =NA,bty="n")
dev.off()

hn<-sce[,sce$patient=="HN1"]
hn$ann1<-as.character(hn$ann1)
hn$ann<-as.character(hn$ann)
load("D:/Bioinfrolf/data/tmpRdata/05.Mono_Macr-tmp.Rdata")
sce<-Mono_Macr
sce$ann<-as.character(sce$ann)
Normal<-sce[,sce$tissue_sub=="Normal"]
Normal<-Normal[,setdiff(colnames(Normal),colnames(hn))]
TumorCore<-sce[,sce$tissue_sub=="TumorCore"]
TumorEdge<-sce[,sce$tissue_sub=="TumorEdge"]

sce<-hn
table(sce$ann1)

#########################scRNAseq Similarity heatmap##############
col_fun = colorRamp2(c(0,1), c("white","#5BBCD6"))#colorRamp2(c(0,1), c("white", "#B40F20"))
jpeg(file="Similarity-nc.jpeg",width=3,height=3,units = "in", res = 1000)
Heatmap(
  b,
  col = col_fun,
  name = 'Similarity',
  show_row_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = F,
  row_title_gp = gpar(fontsize = 16),
  column_title_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16)
)
dev.off()
########################Plot for cellchat######################
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(Seurat)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(SeuratObject)
library(wesanderson)
library(ggsci)
library("RColorBrewer")
col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6","#E2D200","#B40F20","#273046","#FD6467")
col11<-c(col8,c("#7294D4","#E6A0C4","#F2AD00"))
col18<-c(col11c("#DD8D29","#354823","#E58601","#FAD510","#CB2314","#1E1E1E"))

barplot(1:18, col = col18)
load("D:/Bioinfrolf/data/tmpRdata/cpdb-cancer.Rdata")
levels(immune.combined)
table(immune.combined$anno,immune.combined$seurat_clusters)

n<-colnames(normal)
p<-colnames(para)
c<-colnames(core)

immune.combined$tissue<-ifelse(colnames(immune.combined)%in%c,"Core Tumor",ifelse(
  colnames(immune.combined)%in%p,"Peripheral Tumor","Adjacent Normal"
))
table(immune.combined$tissue)

setwd("D:/HCC-SC/cancer/figs")
sce<-immune.combined
table(sce$anno)
table(sce$cluster)
head(sce@meta.data)
umap = sce@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = Idents(sce))

jpeg(file = "umap-all-1.jpg", width = 6, height = 4, units = "in", res = 2000)
DimPlot(sce, reduction = "umap",group.by = "seurat_clusters",label = TRUE, pt.size = .1,cols = col18)
dev.off()

umap = sce@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = sce@meta.data$cluster)

jpeg(file = "umap-all-cluster.jpg", width = 6, height = 4, units = "in", res = 2000)
#pdf(file="umap-all-anno.pdf",width=6,height=4)
ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
  geom_point(size = 0.01, alpha = 1)+
  labs(title ="HCC immune cells",
       x = 'UMAP_1',
       y = "UMAP_2")+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 15,face = 'bold'),
        axis.text = element_text(size=12), # Font size of axis labels.
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
        #legend.position = c(.92, .9),                     # Position of the legend.
        #legend.margin = margin(6, 6, 6, 6),                # Margin of the legend.
        #legend.background = element_rect(size = 0.2, colour = 1)
  )+
  scale_color_manual(values=col11[4:7])
dev.off()

umap = sce@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = sce@meta.data$tissue)

jpeg(file = "umap-all-tissue.jpg", width = 6, height = 4, units = "in", res = 2000)
#pdf(file="umap-all-anno.pdf",width=6,height=4)
ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
  geom_point(size = 0.01, alpha = 1)+
  labs(title ="HCC immune cells",
       x = 'UMAP_1',
       y = "UMAP_2")+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        title = element_text(size = 15,face = 'bold'),
        axis.text = element_text(size=12), # Font size of axis labels.
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_blank(), 
        legend.key.size=unit(0.2, "inches")
        #legend.position = c(.92, .9),                     # Position of the legend.
        #legend.margin = margin(6, 6, 6, 6),                # Margin of the legend.
        #legend.background = element_rect(size = 0.2, colour = 1)
  )+
  scale_color_manual(values=c("#00A08A","#FF0000","#F98400"))
dev.off()
