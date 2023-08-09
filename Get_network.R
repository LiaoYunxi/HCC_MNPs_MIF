suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="a gene expression matrix", metavar="character"),
    make_option(c("--TF"), type="character", default=NULL, 
                help="The Cancer-related TFs obtained through CistromeCancer", metavar="character"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.01, 
                help="pvalue", metavar="numeric"),
    make_option(c("-c", "--corFilter"), type="numeric", default=0.4, 
                help="Select the interactions with |correlation coefficient| > corFilter", metavar="numeric"),
    make_option(c("--name"), type="character", default=NULL, 
                help="Dataset name", metavar="character"),
    make_option(c("--gene1"), type="character", default=NULL,  
                help="target gene 1", metavar="character"),
    make_option(c("--gene2"), type="character", default=NULL,  
                help="target gene 2", metavar="character"),
    make_option(c("--top"), type="integer", default=50,  
                help="top x highest expressed samples", metavar="integer"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="network pic name", metavar="character")
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  file.input = opt$input
  TF =opt$TF
  gene1 = opt$gene1
  gene2 = opt$gene2
  j = opt$top
  pvalueFilter = opt$pvalue
  corFilter = opt$corFilter
  file.output = opt$monole3_output
 
  data <- read.table(file.input,sep="\t")

  # https://cistrome.org/CistromeCancer/CancerTarget/
  TF_list=read.table(TF,sep="\t",header=F)
  tf=intersect(rownames(data),TF_list$V1)
  datatf=data[tf,]

  # filter DEGs
  filterDEG<-function(res){
    res<-na.omit(res[abs(res$logFC)!=Inf,])
    res<-res[res$pValue<0.05,]
    res<-res[order(res$logFC,decreasing = T),]
    rownames(res)<-res$gene
    logFC_cutof = with(res,mean(abs(res$logFC)) + 2*sd(abs(res$logFC)))
    res$result = as.factor(ifelse(res$pValue < 0.05 & abs(res$logFC) >=logFC_cutof,
                                  ifelse(res$logFC >= logFC_cutof ,'UP','DOWN'),'NOT'))
    return(res)
  }

  TM<-read.csv(paste0(name,'-',gene1,"-",j,"-DEG.csv"))
  TS<-read.csv(paste0(name,'-',gene2,"-",j,"-DEG.csv"))
  TM<-filterDEG(TM)
  TS<-filterDEG(TS)

  # shared DEGs
  tg<-intersect(TM$gene[TM$result!="NOT"],TS$gene[TS$result!="NOT"])
  datadg=data[tg,]

  # Select the interactions with |correlation coefficient| > 0.4
  outTab=data.frame()
  for(i in row.names(datatf)){
    if(sd(datatf[i,])>1){
      for(j in row.names(datadg)){
        x=as.numeric(datatf[i,])
        y=as.numeric(datadg[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(TF=i,datadg=j,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(TF=i,datadg=j,cor,pvalue,Regulation="negative"))
        }}}}
  write.table(file=paste0(name,'-',gene1,'-',gene2,"-corResults.txt"),outTab,sep="\t",quote=F,row.names=F) #correlation table

  # Calculate the degree of centrality 
  edges <-outTab[,1:2]
  colnames(edges) <- c("from", "to")
  nodes <- data.frame(name = unique(union(edges$from, edges$to)))
  nodes$type=ifelse(nodes$name%in% TF$V1,"TF","target")

  # Map the gene regulatory networks
  suppressMessages(library(tidygraph,quietly = TRUE))
  suppressMessages(library(ggraph,quietly = TRUE))
 
  Tgraph <- as_tbl_graph(outTab, directed = F,edges =outTab$cor)
  Tg<-Tgraph %>% mutate(centrality = centrality_degree(weights = cor))
 
  pdf(file=paste0(file.output,".pdf"),width =7,height =5)
  ggraph(Tg,layout = 'linear', circular = TRUE) + 
    geom_edge_link(aes(edge_width=cor)) +
    geom_node_point(aes(size = centrality,fill=factor(nodes$type),
                        colour=factor(nodes$type))) +
    geom_node_text(aes(filter= centrality>3,label = name),repel = TRUE)+
    theme_graph()
  dev.off()
 }

main()