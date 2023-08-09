suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Seurat Object of single-cell RNA expression", metavar="character"),
    make_option(c("--x1"), type="numeric", default=4.5, 
                help="root cell x start location", metavar="numeric"),
    make_option(c("--x2"), type="numeric", default=5,  
                help="root cell x end location", metavar="numeric"),
    make_option(c("--y1"), type="numeric", default=-2.25,  
                help="root cell y start location", metavar="numeric"),
    make_option(c("--y2"), type="numeric", default=-1.75,  
                help="root cell y end location", metavar="numeric"),
    make_option(c("-o", "--monole3_output"), type="character", default='monole3_output.Rdata', 
                help="Monole2 Output", metavar="character"),
    make_option(c("-g", "--gene_output"), type="character", default='monocle3.gene.txt', 
                help="Driver Gene Output", metavar="character") 
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if(is.null(opt$input)){
    print_help(opt_parser)
    stop("Seurat Object file must be specified", call.=FALSE)
  }
  
  file.input = opt$input
  x1 = opt$x1
  x2 = opt$x2
  y1 = opt$y1
  y2 = opt$y2
  file.output = opt$monole3_output
  gene.output = opt$gene_output
 
  suppressMessages(library(Seurat))
  suppressMessages(library(monocle))
  suppressMessages(library(tidyverse))

  print('Monole2 Run ...')
  scRNAdata <- readRDS(file.input)
  monocle.matrix=as.matrix(scRNAdata@assays$RNA@data)
  monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
  monocle.sample=as.matrix(scRNAdata@meta.data)
  monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
  monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
  monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
  clusterAnn=as.character(monocle.geneAnn[,2])
  names(clusterAnn)=paste0("cluster",monocle.geneAnn[,1])
 
  scRNAdata.markers <- FindAllMarkers(object = scRNAdata,
                                      only.pos = FALSE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.5)
  sig.markers=scRNAdata.markers[(abs(as.numeric(as.vector(scRNAdata.markers$avg_log2FC)))>0.5 & 
                                 as.numeric(as.vector(scRNAdata.markers$p_val_adj))<0.05),]
 
  data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
  pd<-new("AnnotatedDataFrame", data = monocle.sample)
  fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
  cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
  names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
  pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
  pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)
 
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- setOrderingFilter(cds, scRNAdata.markers$gene)
  cds <- reduceDimension(cds, max_components = 3,verbose = F,reduction_method = 'tSNE')
 
  ## Replace dimensionality reduction data in cds object as the codes of the 18th step.
  cds.embed <- cds@int_colData$reducedDims$TSNE
  int.embed <- Embeddings(scRNAdata, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
 
  cds <- orderCells(cds)
 
  ### for driver gene
  GetBeamGene <- function(branch_point,cds){
    beam_genes <- BEAM(cds, branch_point, verbose = F, cores = detectCores() - 2)
    return(rownames(beam_genes[beam_genes$use_for_ordering=="TRUE"]))
  }
  ## Get beam genes of each branch
  beam_genes1 <- GetBeamGene(1,cds)
  beam_genes2 <- GetBeamGene(2,cds)
  beam_genes3 <- GetBeamGene(3,cds)
  beam_genes<-Reduce(intersect,list(beam_genes1,beam_genes2,beam_genes3))

  ## Get genes associated with differentiation
  diff_test_res <- differentialGeneTest(cds[marker$gene, ], 
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                        cores = detectCores() - 2)
  sig_genes <- row.names(subset(diff_test_res, qval < 0.05))

  gene<-union(beam_genes,sig_genes)

  ### save results
  save(cds, file = file.output)
  write.table(gene,file = gene.output,row.names = F,quote = F)
 }
main()