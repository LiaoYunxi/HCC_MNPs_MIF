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
    make_option(c("-m", "--morans_I"), type="numeric", default=0.1,  
                help="Morans_I", metavar="numeric"),
    make_option(c("-o", "--monole3_output"), type="character", default='monole3_output.Rdata', 
                help="Monole3 Output", metavar="character"),
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
  m = opt$morans_I
  file.output = opt$monole3_output
  gene.output = opt$gene_output
 
  suppressMessages(library(Seurat))
  suppressMessages(library(monocle3))

  print('Monole3 Run ...')
  scRNAdata <- readRDS(file.input)
  data <- GetAssayData(scRNAdata, assay = 'RNA', slot = 'counts')
  cell_metadata <- scRNAdata@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = 50)

  ## Replace dimensionality reduction data in cds object
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(scRNAdata, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)

  ## set root cells
  embed <- data.frame(Embeddings(scRNAdata, reduction = "umap"))
  embed <- subset(embed, UMAP_1 > x1 & UMAP_1 < x2 & UMAP_2 > y1 & UMAP_2 < y2)
  root.cell <- rownames(embed)

  cds <- order_cells(cds, root_cells = root.cell)

  ### Get driver gene
  Track_genesForT<-Track_genes[Track_genes$status=="OK",]
  Track_genesForT<-Track_genesForT[Track_genesForT$p_value<0.05,]
  Track_genesForT<-Track_genesForT[order(Track_genesForT$morans_I,decreasing = T),]

  ### save results
  save(cds, file = file.output)
  write.table(Track_genesForT[Track_genesForT$morans_I>0.1,], file = gene.output,row.names = F,quote = F)
 }

main()