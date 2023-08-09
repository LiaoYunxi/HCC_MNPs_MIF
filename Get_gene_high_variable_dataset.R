suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="a gene expression matrix", metavar="character"),
    make_option(c("--name"), type="character", default=NULL, 
                help="Dataset name", metavar="character"),
    make_option(c("--gene"), type="character", default=NULL,  
                help="target gene", metavar="character"),
    make_option(c("--top"), type="integer", default=50,  
                help="top x highest expressed samples", metavar="integer")
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  file.input = opt$input
  name = opt$name
  gene = opt$gene
  j = opt$top
 
  data <- read.table(file.input,sep="\t")
  treatNum=conNum=j
  grade=c(rep(1,j),rep(2,j))
 
  # Rank patient samples by target gene expression levels 
  OD<-function(data,gene){
    data<-t(data)%>%as.data.frame()
    index<-which(colnames(data)==gene)
    data<-data[order(data[,index],decreasing = F),]
    data<-t(data)%>%as.data.frame()
  }
  rt<-OD(data,gene) # ex.TCGA.tumor is a gene expression matrix

  # Constructed high-variable datasets from the top 50 and tail 50 samples
  data<-rt[,c(1:j,(ncol(rt)-j+1):ncol(rt))]

  # Get DEGs
  suppressMessages(library(optparse,quietly = stats))
  outTab=data.frame()
  for(i in row.names(data)){
      tmp=rbind(expression=data[i,],grade=grade)
      tmp=as.matrix(t(tmp))
      wilcoxTest<-wilcox.test(expression ~ grade, data=tmp)
      conGeneMeans=mean(as.numeric(data[i,1:conNum]))
      treatGeneMeans=mean(as.numeric(data[i,(conNum+1):ncol(data)]))
      logFC=log2(treatGeneMeans)-log2(conGeneMeans)
      pvalue=wilcoxTest$p.value
      outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,
                                treatMean=treatGeneMeans,
                                logFC=logFC,pValue=pvalue))}
  pValue=outTab[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
  outTab=cbind(outTab,fdr=fdr)
  outTab$pValue<-as.numeric(outTab$pValue)
  outTab=outTab[order(outTab$pValue,decreasing = F),]
  write.csv(outTab,file = paste0(name,'-',gene,"-",j,"-DEG.csv"),quote = F)
 }

main()