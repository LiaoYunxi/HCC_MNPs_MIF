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
