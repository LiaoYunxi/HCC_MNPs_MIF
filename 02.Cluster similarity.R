rm(list=ls())
library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)       

load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/05.Mono_Macr_noCCL5.Rdata")
load("D:/Bioinfrolf/Bioinfrolf/tmpRdata/hnsce.Rdata")

glm.predict <-
  function(train.data, train.group, downsample = FALSE, sample.cells = 0, genes.used = NA, test.data, test.group, alpha = 0.99, nfolds = 10) {
    ## Calculate the similarities of the train data and test data.
    ##
    ## Args:
    #' @train.data: A train data matrix with each cell in a column and each gene
    #' in a row.
    #' @train.group: A vector with the same length as the column of train.data.
    #' @downsample: Whether to sample cells in each cluster to the minimum cluster size.
    #' @sample.cells: Sample cells in each group of cells in train data, if 0 do not
    #' sample cells.
    #' @genes.used: Use a subset of genes in both the train and test data.
    #' @test.data: A test data matrix with each cell in a column and each gene
    #' in a row.
    #' @test.group: A vector with the same length as the column of train.data.
    #' @alpha: The elasticnet mixing parameter, with 0≤α≤1, passed to cv.glmnet.
    #' @nfolds: Number of folds, passed to cv.glmnet.
    ##
    ## Returns:
    ## The probability of each cell in the test.data to be predicted as each group.
    require(glmnet)
    require(ComplexHeatmap)
    glm.fits <- list()
    glm.predict <- list()
    if (length(genes.used) > 1) {
      train.data <- train.data[genes.used,]
      test.data <- test.data[genes.used,]
      if (length(genes.used) <= 50) {
        cat("There were less than 50 features used in the training data!\n")
      }
    }
    if (sample.cells == 0 & downsample) {
      sample.cells <- max(50, min(table(train.group)))
    }
    if (sample.cells > 0) {
      ngroup <- length(unique(train.group))
      if (ncol(train.data) >= sample.cells * ngroup) {
        cells_used <- c()
        for (groupi in sort(unique(train.group))) {
          if (length(which(train.group == groupi)) > sample.cells) {
            cells_used <-
              c(cells_used, sample(which(train.group == groupi), sample.cells))
          } else{
            cells_used <- c(cells_used, which(train.group == groupi))
          }
        }
        train.data <- train.data[, cells_used]
        train.group <- train.group[cells_used]
      }
    }
    for (groupi in sort(unique(train.group))) {
      fac <-  factor(train.group == groupi)
      glm.fits[[groupi]] <-
        cv.glmnet(x = t(train.data), fac, offset = getPopulationOffset(fac), 
                  family = 'binomial', intercept = FALSE, 
                  alpha = alpha, nfolds = nfolds, type.measure = 'class'
        )
      glm.predict[[groupi]] <-
        predict(
          object = glm.fits[[groupi]],
          newx = t(test.data),
          newoffset = rep(0, ncol(test.data)),
          s = 'lambda.min'
        )
    }
    glm.predict.df <- data.frame(do.call(cbind, glm.predict))
    colnames(glm.predict.df) <- sort(unique(train.group))
    glm.predict.df.prob <- (1 + exp(-glm.predict.df)) ** -1
    glm.cluster <-
      colnames(glm.predict.df.prob)[apply(glm.predict.df.prob, 1, which.max)]
    glm.predict.mean <-
      apply(glm.predict.df, 2, function(e)
        sapply(split(e, test.group), mean))
    glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
    heatmap <- Heatmap(
      t(glm.predict.mean.prob),
      name = 'Predicted\nSimilarity',
      column_title = 'test data',
      row_title = 'train data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
    return(
      list(
        test.group = test.group,
        logits = glm.predict.df,
        probability = glm.predict.df.prob,
        cluster = glm.cluster,
        heatmap = heatmap
      )
    )
  }

getPopulationOffset = function(y) {
  ## Calculate the offset value used in glm.predict.
  if (!is.factor(y))
    y = factor(y)
  if (length(levels(y)) != 2)
    stop("y must be a two-level factor")
  off = sum(y == levels(y)[2]) / length(y)
  off = log(off / (1 - off))
  return(rep(off, length(y)))
}

Mono_Macr$ann<-Mono_Macr$ann[drop=T]
old.cluster.name <- levels(hn$ann1)
new.cluster.name <- c("KC-VCAM1","MDM-HLA-DQA1","MDM-APOE","MDM-THBS1","MDM-CD4","MDM-VCAN")
hn$ann1 <- plyr::mapvalues(x = hn$ann1, from = old.cluster.name, to = new.cluster.name)
hn$ann1<-factor(hn$ann1,levels = c("KC-VCAM1","MDM-APOE","MDM-HLA-DQA1","MDM-CD4","MDM-THBS1","MDM-VCAN"))
ada1 <- Mono_Macr
ada2 <- hn

ada1_df <- ada1$ann
ada2_df <- ada2$ann1

gene1 <- rownames(ada1@assays$RNA@data)
gene2 <- rownames(ada2@assays$RNA@data)

ada1_exp <- ada1@assays$RNA@data
ada1_exp[1:5,1:5]
#rownames(ada1_exp) <- rownames(py_to_r(ada1$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- ada2@assays$RNA@data
ada2_exp[1:5,1:5]
#rownames(ada2_exp) <- rownames(py_to_r(ada2$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2))

res <- glm.predict(ada1_exp, ada1_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada2_exp, ada2_df, alpha = 0.99, nfolds = 10)

glm.predict.mean <-
  apply(res$logits, 2, function(e)
    sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
glm.predict.mean.prob<-t(glm.predict.mean.prob)
glm.predict.mean.prob<-glm.predict.mean.prob[c(2,1,4,3,5,7,6),]
library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("light grey","white", "#FF0000"))#c("#e9e9e9","white", "red")

jpeg(file="clusters-Similarity-1.jpeg",width=8,height=6,units = "in", res = 1000)
Heatmap(
  glm.predict.mean.prob,
  col = col_fun,
  name = 'Predicted\nSimilarity',
  column_title = 'Healthy Liver MNPs',
  row_title = 'HCC MNPs',
  show_row_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title_gp = gpar(fontsize = 16),
  column_title_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16)
)
dev.off()
