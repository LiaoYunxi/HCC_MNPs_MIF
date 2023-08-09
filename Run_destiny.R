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
    make_option(c("-o", "--output"), type="character", default='destiny_output.Rdata', 
                help="destiny Output", metavar="character"),
    make_option(c("-g", "--gene_output"), type="character", default='destiny.gene.txt', 
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
  file.output =opt$output
  gene.output = opt$gene_output
 
  suppressMessages(library(Seurat))
  suppressMessages(library(destiny))

  print('destiny Run ...')
 
  scRNAdata <- readRDS(file.input)
  matrix=as.matrix(scRNAdata@assays$RNA@data)
  pData<-scRNAdata@meta.data
  adf<-as(pData, "AnnotatedDataFrame")
  myExpressionSet <- ExpressionSet(assayData=matrix,
                                   phenoData=adf,
                                   annotation="hgu95av2")
  sigmas <- find_sigmas(myExpressionSet, verbose = FALSE,
                        censor_val = 10, censor_range = c(10, 40))
  dm <- DiffusionMap(myExpressionSet, sigmas, verbose = FALSE,
                       censor_val = 10, censor_range = c(10, 40))
  ## set root cells
  embed <- data.frame(Embeddings(scRNAdata, reduction = "umap"))
  embed <- subset(embed, UMAP_1 > x1 & UMAP_1 < x2 & UMAP_2 > y1 & UMAP_2 < y2)
  root.cell <- rownames(embed)

  DPT<-DPT(dm, tips = which(rownames(pData)%in%sample(root.cell,3)))
  myExpressionSet$pseudotime_diffusionmap<-DPT$dpt #rank(DPT$dpt)
 
  # dm is a DiffsionMap object after running destiny pipline
  gr <- gene_relevance(dm)
  gms <-plot(gr, iter_smooth = 0)
  save(gms,gr,myExpressionSet,dm, file = file.output)
  write.table(names(gms$scores),file = gene.output,quote = F,row.names = F)
 }

main()
