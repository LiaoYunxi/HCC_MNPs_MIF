rm(list = ls())
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(wesanderson)
library(ggsci)
library(RColorBrewer)
library(patchwork)
options(stringsAsFactors = FALSE)

if(color){
  col3<-wes_palettes$Darjeeling1[c(1,4,2)]
  col4<-wes_palettes$Darjeeling1[c(1,4,2,5)]
  col6<-c(wes_palettes$Darjeeling1[c(1,2,4,5)],
          wes_palettes$FantasticFox1[c(2,5)])
  annocol1=c("Macr-APOE" ="#FF0000","Macr-APOC1"="#F98400","Mono-VCAN"="#273046",
          "Mono-FCGR3A"="#FD6467","Macr-MIF10"="#E2D200","Macr-VCAM1"="#00A08A","Macr-SPP1"="#5BBCD6",
          "Hepatocytes"="#E6A0C4","Fibroblasts"="#F2AD00",'Endothelial Cells'="#7294D4")

  col8<-c("#FF0000","#00A08A","#F98400","#5BBCD6","#E2D200","#273046","#FD6467")
  
  col11<-c(col8,c("#7294D4","#E6A0C4","#F2AD00"))

  col14<-c(col8,c("#354823","#F2AD00","#E6A0C4","#C6CDF7","#D8A499","#7294D4"))
  barplot(1:10,col = col11)
}
if(cell combinding){
  groupSize <- as.numeric(table(cellchat@idents))
  vertex.receiver = seq(1,7) #a numeric vector
  
  load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/cpdb-cancer.Rdata")
  table(immune.combined$anno)
  immune.combined<-subset(x = immune.combined,subset = anno!="Macr-MIF5")
  immune.combined$anno<-immune.combined$anno[drop=T]
  
  normal<-subset(x = normal,subset = anno!="Macr-MIF5")
  normal$anno<-normal$anno[drop=T]
  table(normal$anno)
  normal<-subset(x = normal,subset = anno!="Macr-MIF5")
  normal$anno<-normal$anno[drop=T]
  table(normal$anno)
  normal<-subset(x = normal,subset = anno!="Macr-MIF5")
  normal$anno<-normal$anno[drop=T]
  table(normal$anno)
  
  save(list = ls(),file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/cpdb-cancer-noMIF5.Rdata")
  rm(immune.combined)
  
}
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/cpdb-cancer-noMIF5.Rdata")
data.input = normal@assays$RNA@data #normalized data matrix
meta = normal@meta.data # a dataframe with rownames containing cell mata data
#cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$anno) # check the cell labels
head(normal@meta.data)

identity = data.frame(group = normal$anno, row.names = names(normal$anno)) 
unique(identity$group) 

cellchat <- createCellChat(object = data.input,meta =normal@meta.data,group.by = "anno")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "anno")
levels(cellchat@idents)
#str(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)    # Show the structure of the database
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost

future::plan("multinormal", workers = 12) # do normalllel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

setwd("D:/Bioinfrolf/HCC-SC/cancer/cellchat.21.12/normal")
#df.net <- subsetCommunication(cellchat)
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.net <- subsetCommunication(cellchat, sources.use = c(8:10), targets.use = c(1:7)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat@LR$LRsig$pathway_name%>%table()
cellchat@netP$pathways
head(cellchat@LR$LRsig)
vertex.receiver = seq(1,7) #a numeric vector

##########################circle#########################
dev.off()
library(patchwork)
jpeg(file = "cellchat-net-normal-SS.jpg", width = 10, height = 6, units = "in", res = 2000)
par(mfrow = c(1,2), xpd=TRUE)
p1<-netVisual_circle(cellchat@net$count, vertex.weight = groupSize,color.use=col11, weight.scale = T, label.edge= F, title.name = "Number of interactions")
p2<-netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,color.use=col11, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
p1+p2
dev.off()

mat <- cellchat@net$weight
jpeg(file = "cellchat-shell-normal-SS.jpg", width = 16, height = 12, units = "in", res = 1000)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,color.use=col11,weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.use = col11,color.heatmap ="YlGnBu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.use = col11,color.heatmap ="YlGnBu")
jpeg(file = "cellchat-signalingRole-hp-normal-SS.jpg", width = 10, height = 5, units = "in", res = 1000)
ht1 + ht2
dev.off()
?netAnalysis_signalingRole_heatmap()
######################################chord plot#################################
setwd("D:/Bioinfrolf/HCC-SC/cancer/cellchat.21.12/normal")
jpeg(file = "cellchat-bubble-normal-Hepatocytes-SS.jpg", width = 5, height = 4, units = "in", res = 1000)
netVisual_bubble(cellchat, sources.use = 9,
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)
dev.off()

jpeg(file = "cellchat-bubble-normal-Fibroblast-SS.jpg", width = 5, height = 4, units = "in", res = 1000)
netVisual_bubble(cellchat, sources.use = 10,
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)
dev.off()

jpeg(file = "cellchat-bubble-normal-Endothelial Cells-SS.jpg", width = 5, height = 4, units = "in", res = 1000)
netVisual_bubble(cellchat, sources.use = 11,
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)
dev.off()

jpeg(file = "cellchat-chord-normal-Endothelial Cells-SS.jpg", width = 5, height = 5, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = 11, targets.use = c(1:7), color.use=col11,show.legend = F,lab.cex = 1.5)
dev.off()

jpeg(file = "cellchat-chord-normal-Hepatocytes-SS.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = 9, targets.use = c(1:7), lab.cex = 1.5, color.use=col11,show.legend = F)
dev.off()

jpeg(file = "cellchat-chord-normal-Fibroblast-SS.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = 10, targets.use = c(1:7), lab.cex = 1, color.use=col11,show.legend = F)
dev.off()

jpeg(file = "cellchat-chord-netP-normal-SS.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = c(8:10), targets.use = c(1:7), slot.name = "netP", lab.cex = 1,legend.pos.x = 10,color.use=col11)
dev.off()

jpeg(file = "cellchat-chord-normal-SS.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = c(8:10), targets.use = c(1:7),lab.cex = 1,legend.pos.x = 10,color.use=col11)
dev.off()

jpeg(file = "cellchat-chord-normal-SS-reverse.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = c(1:7), targets.use = c(8:10),lab.cex = 1,legend.pos.x = 10,color.use=col11)
dev.off()
#netVisual_chord_gene(cellchat, sources.use = c(8:10), targets.use = c(1:7), signaling = c("MIF","MIF"),legend.pos.x = 8)

#######################################NMF############
library(NMF)
library(ggalluvial)
selectK(cellchat,slot.name = "netP",pattern = "incoming")
dev.off()
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

jpeg(file = "cellchat-netAnalysis_river-normal-in-SS.jpg", width = 6, height = 4, units = "in", res = 2000)
netAnalysis_river(cellchat, pattern = "incoming",color.use=col11)
dev.off()
jpeg(file = "cellchat-netAnalysis_dotr-normal-in-SS.jpg", width = 5, height = 4, units = "in", res = 2000)
netAnalysis_dot(cellchat, pattern = "incoming",color.use=annocol1)
dev.off()

selectK(cellchat,slot.name = "netP",pattern = "outgoing")
nPatterns = 4
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)  
dev.off()

jpeg(file = "cellchat-netAnalysis_river-normal-out-SS-5.jpg", width = 6, height = 4, units = "in", res = 2000)
netAnalysis_river(cellchat, pattern = "outgoing",color.use=col11)
dev.off()
jpeg(file = "cellchat-netAnalysis_dotr-normal-out-SS.jpg", width = 5, height = 4, units = "in", res = 2000)
netAnalysis_dot(cellchat, pattern = "outgoing",color.use=annocol1)
dev.off()



####################################################################
noNAnetCluster<-function (object, slot.name = "netP", type = c("functional", 
                                                               "structural"), comparison = NULL, k = NULL, methods = "kmeans", 
                          do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, 
                          k.eigen = NULL) {
  type <- match.arg(type)
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Classification learning of the signaling networks for a single dataset", 
        "\n")
  }
  else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Classification learning of the signaling networks for datasets", 
        as.character(comparison), "\n")
  }
  comparison.name <- paste(comparison, collapse = "-")
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  data.use <- na.omit(Y)
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use, k, nstart = 10)$cluster
    }
    else {
      N <- nrow(data.use)
      kRange <- seq(2, min(N - 1, 10), by = 1)
      if (do.parallel) {
        future::plan("multiprocess", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(test = future::nbrOfWorkers() == 
                            1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
      results = my.sapply(X = 1:length(kRange), FUN = function(x) {
        idents <- kmeans(data.use, kRange[x], nstart = 10)$cluster
        clusIndex <- idents
        adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, 
                                                   clusIndex, FUN = "==")), nrow = N, ncol = N)
        return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
      }, simplify = FALSE)
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce("+", adjMat)/length(kRange)
      res <- computeEigengap(as.matrix(CM))
      numCluster <- res$upper_bound
      clusters = kmeans(data.use, numCluster, nstart = 10)$cluster
      if (do.plot) {
        gg <- res$gg.obj
        ggsave(filename = paste0("estimationNumCluster_", 
                                 fig.id, "_", type, "_dataset_", 
                                 comparison.name, ".pdf"), plot = gg, 
               width = 3.5, height = 3, units = "in", 
               dpi = 300)
      }
    }
  }
  else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D) - A
    L <- diag(D^-0.5) %*% L %*% diag(D^-0.5)
    evL <- eigen(L, symmetric = TRUE)
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[, (ncol(evL$vectors) - k.eigen + 1):ncol(evL$vectors)]
    clusters = kmeans(Z, k, nstart = 20)$cluster
  }
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
    methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
  return(object)
}

load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/cellchat-normal-SS.Rdata")
normal<-cellchat
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/cellchat-para-SS.Rdata")
para<-cellchat
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/cellchat-core-SS.Rdata")
core<-cellchat

cellchat<-para

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
cellchat <-noNAnetCluster(cellchat, type = "functional")

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <-noNAnetCluster(cellchat, type = "structural")

object<-cellchat
slot.name = "netP"
type = "structural"#"functional"#
color.use = NULL
pathway.remove = NULL
pathway.remove.show = TRUE
dot.size = c(1, 5)
label.size = 3.5
dot.alpha = 0.5
xlabel = "Dim 1"
ylabel = "Dim 2"
title = NULL
font.size = 10
font.size.title = 12
do.label = T
show.legend = T
show.axes = T
comparison <- "single"
comparison.name <- paste(comparison, collapse = "-")

Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
Groups <- methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]]
prob <- methods::slot(object, slot.name)$prob
similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
pathway.remove <- rownames(similarity)[which(colSums(similarity) ==1)]
pathway.remove.idx <- which(dimnames(prob)[[3]] %in% pathway.remove)
#prob <- prob[, , -pathway.remove.idx]
prob_sum <- apply(prob, 3, sum)
df <- data.frame(x = Y[, 1], y = Y[, 2], Commun.Prob. = prob_sum/max(prob_sum), 
                 labels = as.character(unlist(dimnames(prob)[3])))
df$Groups<-NA
df$Groups[which(rownames(df)%in%names(Groups))]<-Groups
df$Groups<-as.factor(df$Groups)
col8<-c( "#E2D200","#5BBCD6","#B40F20","#FD6467","#F98400","#273046")
color.use <- col8[1:length(unique(Groups))]

gg <- ggplot(data = na.omit(df), aes(x, y)) + geom_point(aes(size = Commun.Prob.,fill= Groups,color=Groups), shape = 21) + 
  CellChat_theme_opts() + 
  theme(text = element_text(size = font.size), legend.key.height = grid::unit(0.15,"in")) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  labs( x = "UMAP1", y = "UMAP2") + 
  theme(plot.title = element_text(size = font.size.title, face = "plain")) + 
  scale_size_continuous(limits = c(0, 1), range = dot.size, breaks = c(0.1, 0.5, 0.9)) + 
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5),
        title = element_text(size = 15),
        legend.text =element_text(size=10),  # Font size of legend labels.
        legend.title = element_text(size=12), 
        legend.key.size=unit(0.2, "inches")
  )+theme(
    axis.line = element_line(colour = 'black', size = 0.5), 
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12)
  )
gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use,alpha = dot.alpha), drop = FALSE)
gg <- gg + scale_colour_manual(values = color.use, drop = FALSE)
gg <- gg + ggrepel::geom_label_repel(mapping = aes(label = labels,colour = Groups), size = label.size, show.legend = F, 
                                    segment.size = 0.2, segment.alpha = 0.5)
# gg <- gg + annotate(geom = "text", 
#                     label = paste("Isolate pathways: ",paste(c(rownames(df)[!rownames(df)%in%names(Groups)],pathway.remove), collapse = ", ")), 
#                     x = -Inf,y = Inf, hjust = 0, vjust = 1, 
#                     size = label.size,fontface = "italic")


jpeg(file = "cellchat-s_Similarity-para.jpg", width = 4, height =2.5, units = "in", res = 1000)
#netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
gg
dev.off()


save(cellchat,col11,annocol1,normal,file = "D:/Bioinfrolf/Bioinfrolf/tmpRdata/cellchat-normal-SS.Rdata")

#######################pathway plot################################
setwd("D:/Bioinfrolf/HCC-SC/cancer/cellchat.21.12/normal/MIF")
pathways.show <- c("MIF")
# ?netVisual_aggregate()
# ?netAnalysis_contribution()

# Hierarchy plot
jpeg(file = "cellchat-MIF-Hierarchy-normal-SS.jpg", width = 10, height = 5, units = "in", res = 2000)
netVisual_aggregate(cellchat, signaling= pathways.show,vertex.receiver = vertex.receiver,layout = "hierarchy",
                    vertex.weight = groupSSSSSSSize,color.use=col11,
                    show.legend = TRUE,
                    pt.title = 15,title.space=2,
                    vertex.label.cex = 0.5,vertex.size.max = 20,
                    small.gap=1)
dev.off()

# Circle plot
jpeg(file = "cellchat-MIF-Circle-normal-SS.jpg", width = 6, height = 6, units = "in", res = 2000)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", 
                    vertex.size = groupSize,color.use=col11,
                    show.legend = F,
                    pt.title = 15,title.space=5,
                    vertex.label.cex = 1,vertex.size.max = 20,
                    small.gap=1
)
dev.off()

jpeg(file = "cellchat-MIF-contribution-normal-SS.jpg", width = 5, height = 4, units = "in", res = 2000)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

group.cellType <- c(rep("MNPs", 7),rep("others", 3)) 
names(group.cellType) <- levels(cellchat@idents)
jpeg(file = "cellchat-MIF-chord_cell-normal-SS1.jpg", width = 6, height = 6, units = "in", res = 2000)
par(mfrow=c(1,1))
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"),
                     color.use=col11,
                     show.legend = F,
                     small.gap=1)
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",
#                     vertex.size = groupSize,color.use=col11,
#                     show.legend = F,
#                     pt.title = 15,title.space=5,
#                     vertex.label.cex = 1,vertex.size.max = 20,
#                     small.gap=1
#                     )
dev.off()

# ??netAnalysis_signalingRole_scatter()
# ?netAnalysis_computeCentrality()


jpeg(file = "cellchat-MIF-signalingRole_network-normal-SS.jpg", width = 8, height = 4, units = "in", res = 2000)
netAnalysis_signalingRole_network(cellchat, slot.name = "netP",signaling=pathways.show, 
                                  #color.heatmap = "YlOrRd",
                                  color.use=col11)
dev.off()
jpeg(file = "cellchat-MIF-signalingRole_heatmap-normal-SS.jpg", width = 6, height = 5, units = "in", res = 2000)
netVisual_heatmap(cellchat, signaling = pathways.show,  color.use=col11,color.heatmap = "YlOrRd")
#netAnalysis_signalingRole_heatmap(cellchat, slot.name = "netP",signaling=pathways.show,color.heatmap = "YlOrRd")
dev.off()

jpeg(file = "cellchat-MIF-signalingRole_scatter-normal-SS.jpg", width = 6, height = 4, units = "in", res = 2000)
netAnalysis_signalingRole_scatter(cellchat, slot.name = "netP",signaling=pathways.show,
                                  color.use=col11)
dev.off()

jpeg(file = "cellchat-chord-normal-MIF-SS.jpg", width = 7, height = 7, units = "in", res = 1000)
netVisual_chord_gene(cellchat, sources.use = c(8:10), targets.use = c(1:7), signaling = c("MIF"),lab.cex = 1,legend.pos.x = 10,color.use=col11)
dev.off()

jpeg(file = "cellchat-GeneExpression-normal-MIF-SS.jpg", width = 6, height = 4, units = "in", res = 1000)
plotGeneExpression(cellchat, signaling = "MIF",color.use=col11)
dev.off()

#ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MIF", "MIF"))
jpeg(file = "cellchat-signalingRole-hp-MIF-normal-SS.jpg", width = 5, height = 4, units = "in", res = 1000)
netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MIF"))
dev.off()

gg1 <- netAnalysis_signalingRole_scatter(cellchat,color.use=col11)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MIF"),color.use=col11)
jpeg(file = "cellchat-signalingRole-MIF-normal-SS.jpg", width = 10, height = 5, units = "in", res = 1000)
gg1 + gg2
dev.off()