
suppressPackageStartupMessages({
  library("vcfR")
  library("ape")
  library("optparse")
  library("glue")
  library("cluster")    # clustering algorithms
  library("adegenet")
})

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
  gind.mat <- as.matrix(gind$tab)
  gind.mat[is.na(gind.mat)] <- 0
  gind.mat[gind.mat==2] <- 1
  D <- dist(gind.mat)
  return(D)
}

mp2 <- function(D) {

  tree <- nj(D)
  plot.phylo(tree, type = "unrooted", edge.width = 1, font = 1, lab4ut = "axial", show.tip.label = TRUE, cex = 1)
  write.tree(tree, file = "snp_tree.nwk")

  for ( s in attributes(D)$Labels ) {
    fout <- glue("./{s}.png")
    png(fout, width = 700, height = 700)

    plot.phylo(tree, type = "unrooted", edge.width = 1, font = 1, lab4ut = "axial", show.tip.label = TRUE, cex = 1, tip.color = c("red")[(tree$"tip.label"!=s)+1])

    edgelabels(tex=round(tree$edge.length,1), bg=rgb(.8,.8,1,.8), cex=0.75)
    dev.off()
  }
}

makePCAplot <- function(D) {
  png("./pca.png", width = 700, height = 700)
  pco1 <- dudi.pco(D, scannf=FALSE, nf=2)
  s.label(pco1$li*1.1, clab=0, pch="")
  colorplot(pco1$li, pco1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
  abline(v=0,h=0,col="grey", lty=2)
  title("Principal Coordinate Analysis\n-based on SNP distances-")
  dev.off()
}

kmc <- function(gind) {
  png("kmc.png", width = 700, height = 700)

  ## format data
  gind.mat <- as.matrix(gind$tab)
  gind.mat[is.na(gind.mat)] <- 0
  gind.mat[gind.mat==2] <- 1
  filtered.mat <- gind.mat[ , which(apply(gind.mat, 2, var) != 0)]

  ## automate assigning k
  optimal_k <- fviz_nbclust(filtered.mat, kmeans, method = "silhouette")
  k <- which.max(optimal_k$data$y)

  ## generate and plot
  kmeans_result <- kmeans(filtered.mat, centers = k, nstart = 25)
  p2 <- fviz_cluster(kmeans_result, data = filtered.mat, geom="point", main=glue("k={k}"))

  grid.arrange(p2, nrow = 1)
  dev.off()
}

main <- function(){
  ## Main function
  opt <- parseArgs()
  vcf=opt$vcf
  D <- vcf2dist(vcf)
  mp2(D)

  makePCAplot(D)
  # kmc(gind)

}

main()
