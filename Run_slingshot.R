suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Seurat Object of single-cell RNA expression", metavar="character"),
    make_option(c("-r", "--root_cell"), type="character", default=NULL, 
                help="root cell", metavar="character"),
    make_option(c("--group"), type="character", default=NULL, 
                help="clueter label", metavar="character"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.01, 
                help="pvalue", metavar="numeric"),
    make_option(c("-n", "--number_gene"), type="integer", default=1000, 
                help="the number pf driver genes", metavar="numeric"),
    make_option(c("-g", "--gene_output"), type="character", default='slingshot.gene.txt', 
                help="Driver Gene Output", metavar="character") 
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if(is.null(opt$input)){
    print_help(opt_parser)
    stop("Seurat Object file must be specified", call.=FALSE)
  }
  
  file.input = opt$input
  root_cell = opt$root_cell
  group = opt$group
  pvalue = opt$pvalue
  num = opt$number_gene
  gene.output = opt$gene_output
 
  suppressMessages(library(Seurat))
  suppressMessages(library(slingshot))
  suppressMessages(library(tradeSeq))

  print('slingshot Run ...')
 
  scRNAdata <- readRDS(file.input)
  test.count=scRNAdata@assays$RNA@counts
  test.data=scRNAdata@assays$RNA@data#GetAssayData(scRNAdata,slot = "data")
  test.log=log2(test.count+1)
  clusters<-scRNAdata@meta.data$seurat_clusters
  sce<- SingleCellExperiment(assays = List(counts = test.count,
                                           data=test.data,logcounts=test.log))
  rm(scRNAdata)

  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x > 0) >= 100
  })
  sce <- sce[geneFilter, ]
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  pca <- scRNAdata@reductions$harmony
  rd1 <- pca@cell.embeddings[,1:2]
  rd2<-scRNAdata@reductions$umap@cell.embeddings
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  SingleCellExperiment::reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

  sce@colData$ann<-scRNAdata@meta.data$ann
  sce <- slingshot::slingshot(sce, clusterLabels = group, 
                              reducedDim = 'UMAP',start.clus = root_cell)

  ### Get driver Gene
  # fit negative binomial GAM
  sce <- fitGAM(sce)
  # test for dynamic expression
  ATres <- associationTest(sce)
  genedf <- ATres[ATres$pvalue<pvalue, ]
  write.table(rownames(genedf[1:num,]),
              file = gene.output,quote = F,row.names = F)
 }

main()