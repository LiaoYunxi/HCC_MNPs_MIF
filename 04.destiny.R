library(Biobase)
library(destiny)
library(tidyverse)
library(monocle)
library(ggridges)
library(RColorBrewer)
library(scales)
library(ggbeeswarm)
#BiocManager::install("ouija")
rm(list = ls())
setwd("D:/Bioinfrolf/Bioinfrolf/destiny/")
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr-tmp.Rdata")
Mono_Macr<-subset(x = Mono_Macr,subset = ann!="Macr-CCL5")
Myeloid<-Mono_Macr
data(guo)
str((guo))
guo@phenoData
guo$num_cells
Mono_Macr@meta.data%>%head()
matrix=as.matrix(Myeloid@assays$RNA@data)
pData<-Mono_Macr@meta.data
metadata <-  data.frame(labelDescription=c("tissue_sub",
                                           "ann", 
                                           "tissue"),
                        row.names=c("tissue_sub", "ann", "tissue"))

adf<-as(pData, "AnnotatedDataFrame")
#adf <- new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
myExpressionSet <- ExpressionSet(assayData=matrix,
                                 phenoData=adf,
                                 annotation="hgu95av2")
myExpressionSet$ann
parallel::detectCores() - 2
future::plan("multiprocess", workers = 14) 
?DiffusionMap()
dm_guo <- DiffusionMap(myExpressionSet, verbose = FALSE,
                       #censor_val = 10, censor_range = c(10, 40),
                       n_pcs = 50)
plot(dm_guo,
     col = myExpressionSet$ann, pch = 20)

sigmas <- find_sigmas(myExpressionSet, verbose = FALSE,
                      censor_val = 10, censor_range = c(10, 40))

par(lwd = 3)
plot(sigmas,
     col           = palette()[[1]],
     col_highlight = palette()[[4]],
     col_line      = palette()[[6]])

dm_guo_global <- DiffusionMap(myExpressionSet, sigmas, verbose = T,
                              censor_val = 10, censor_range = c(10, 40),
                              n_pcs = 50)
?plot()

rm(Mono_Macr,matrix)

DPT<-DPT(dm_guo_global, tips = random_root(dm_guo_global))
gr <- gene_relevance(dm_guo_global)

jpeg(file = "destiny_Pseudotime_DPT.jpg", width = 6, height = 4, units = "in", res =1000)
plot.DPT(DPT,
         col_by = "dpt",
         col_path = rev(palette()),
         col_tip = "red",
         col = NULL,
         root=3)
dev.off()

save(list = ls(),file ="destiny-TMP.Rdata" )
dm_guo_global@eigenvectors%>%dim()
dm_guo_global@eigenvectors%>%colnames()

detailed_cell_type_color <- c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
                              "Myeloid-CCL5"="#B40F20","Mono-FCGR3A"="#FD6467","Macr-CXCL10"="#E2D200",
                              "Macr-SPP1"="#5BBCD6","Macr-VCAM1"="#00A08A")

plot(eigenvalues(dm_guo_global), ylim = 0:1, pch = 20,
     xlab = "Diffusion component (DC)", ylab = "Eigenvalue")

jpeg(file = "destiny_DiffusionMap_3d.jpg", width = 6, height = 4, units = "in", res =1000)
plot(dm_guo,
     col = myExpressionSet$ann, pch = 20,
     legend.main = "Cell cluster")
dev.off()

tmp_DPT<-data.frame(DC1=DPT@dm@eigenvectors[,1],
                DC2=dm_guo@eigenvectors[,2],
                Timepoint=myExpressionSet$ann)
tmp<-data.frame(DC1=dm_guo@eigenvectors[,1],
                       DC2=dm_guo@eigenvectors[,2],
                       Timepoint=myExpressionSet$ann)
tmp_global<-data.frame(DC1=dm_guo_global@eigenvectors[,1],
                DC2=dm_guo_global@eigenvectors[,2],
                Timepoint=myExpressionSet$ann)

jpeg(file = "destiny_DiffusionMap.jpg", width = 6, height = 4, units = "in", res =1000)
ggplot(tmp,aes(x=DC1,y=DC2,color=Timepoint))+
        geom_point()+scale_color_manual(values = detailed_cell_type_color)+
        theme_classic()
dev.off()

dd=data.frame(row.names = rownames(dm_guo_global@eigenvectors),pseudotime=dm_guo_global@eigenvectors[,2])
dd_r=data.frame(row.names = rownames(dm_guo_global@eigenvectors),pseudotime=rank(dm_guo_global@eigenvectors[,2]))
myExpressionSet$pseudotime_diffusionmap<-dm_guo_global@eigenvectors[,2]
myExpressionSet$pseudotime_diffusionmap<-rank(dm_guo_global@eigenvectors[,2])

jpeg(file = "destiny_Pseudotime_density.jpg", width = 6, height = 4, units = "in", res =1000)
ggplot(as.data.frame(myExpressionSet@phenoData@data),
       aes(x=pseudotime_diffusionmap,y=ann,fill=ann))+
        geom_density_ridges(scale=1) +
        geom_vline(xintercept = c(5,10),linetype=2)+
        scale_y_discrete("")+
        theme_minimal()+
        theme(
                panel.grid = element_blank()
        )+ scale_fill_manual(values = detailed_cell_type_color)+
        theme(axis.text = element_text(size=12), # Font size of axis labels.
              legend.text =element_text(size=10),  # Font size of legend labels.
              legend.title = element_blank(), 
              legend.key.size=unit(0.2, "inches")
        )
dev.off()
save(dd,dd_r,file = "destinyPseudotime.Rdata")
#library(rgl)
#plot3d(eigenvectors(dm_guo)[, 1:3], col = guo$num_cells)