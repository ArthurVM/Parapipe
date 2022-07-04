library("vcfR")
library("adegenet")
library("ape")
library("optparse")
library("glue")
library("wordcloud")
# library("tidyverse")  # data manipulation
library("cluster")    # clustering algorithms
library("factoextra") # clustering algorithms & visualization
library("gridExtra")

parseArgs <- function() {
  ## Simple argument parser
  
  option_list = list(
    make_option(c("-v", "--vcf"), type="character", default=NULL, 
                help="merged vcf file", metavar="character")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  if ( is.null(opt$vcf) ){
    print_help(opt_parser)
    stop("All arguements need to be supplied.\n", call.=FALSE)
  }
  
  return(opt)
}

vcf2dist <- function(vcf) {
  gind <- vcfR2genind(read.vcfR( vcf , verbose = FALSE ))
  row.names(gind@tab)
  
  return(gind)
}

makePhylo <- function(D) {
  tre <- nj(D)
  ape::write.tree(tre, file="./tree.nwk")
  
  # plot(tre, type="unrooted", edge.w=2)
  # edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
  
  for ( s in tre$"tip.label" ) {
    fout <- glue("./{s}.png")
    png(fout, width = 700, height = 700)
    plot(tre, type="unrooted", edge.w=2, tip.color = c("red")[(tre$"tip.label"!=s)+1])
    edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
    dev.off()
  }
}

makePCAplot <- function(D) {
  pco1 <- dudi.pco(D, scannf=FALSE,nf=2)
  s.label(pco1$li*1.1, clab=0, pch="")
  textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li),
           cex=1.0, new=FALSE, xpd=TRUE)
  title("Principal Coordinate Analysis\n-based on SNP distances-")
}

kmeanscluster <- function(gind) {
  gind.mat <- as.matrix(gind$tab)
  gind.mat[is.na(gind.mat)] <- 0
  gind.mat[gind.mat==2] <- 1
  filtered.mat <- gind.mat[ , which(apply(gind.mat, 2, var) != 0)]
  wss.p <- fviz_nbclust(filtered.mat, kmeans, method = "wss")
  grid.arrange(wss.p, nrow = 1)
  sc.km <- kmeans(filtered.mat, centers = 2, nstart = 5)
  p2 <- fviz_cluster(sc.km, data = filtered.mat, labelsize=0, main="k=2")
  grid.arrange(p2, nrow = 1)
}

main <- function(){
  ## Main function
  opt <- parseArgs()
  # vcf = "/home/amorris/BioInf/Parapipe/Cparvum_AH.vcf"
  gind <- vcf2dist(opt$vcf)
  D <- dist(tab(gind))
  makePhylo(D)
  
  # makePCAplot(D)
  # kmeanscluster(gind)
  
}

main()