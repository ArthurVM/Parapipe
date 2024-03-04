suppressPackageStartupMessages({
  library("vcfR")
  library("ape")
  library("optparse")
  library("glue")
  # library("tidyverse")  # data manipulation
  library("cluster")    # clustering algorithms
  library("adegenet")
})

parseArgs <- function() {
  ## Simple argument parser

  option_list = list(
    make_option(c("-a", "--allele_matrix"), type="character", default=NULL,
                help="allele matrix csv file", metavar="character")
  );

  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);

  if ( is.null(opt$allele_matrix) ){
    print_help(opt_parser)
    stop("All arguements need to be supplied.\n", call.=FALSE)
  }

  return(opt)
}

readAlleleMatrix <- function(allele_matrix) {
  data <- read.csv(allele_matrix, header = TRUE, row.names = 1)
  D <- as.dist(dist(data, method = "euclidean"))
  return(D)
}

mp2 <- function(D) {

  hclust_result <- hclust(D, method = "complete")

  phylo_tree <- as.phylo(hclust_result)
  unrooted_tree <- ladderize(phylo_tree)

  for ( s in unrooted_tree$tip.label ) {
    fout <- glue("./{s}.png")
    png(fout, width = 700, height = 700)

    # plot.phylo(tree, type = "unrooted", edge.width = 1, font = 1, lab4ut = "axial", show.tip.label = TRUE, cex = 1, tip.color = c("red")[(tree$tip.label!=s)+1])
    plot(unrooted_tree, sub = NULL, type = "unrooted", lab4ut = "axial", cex = 0.8, tip.color = c("red")[(unrooted_tree$tip.label!=s)+1])
    nodelabels(pch = 16, col = "blue", cex = 1, srt = 90)

    edgelabels(tex=round(unrooted_tree$edge.length,1), bg=rgb(.8,.8,1,.8), cex=0.75)
    dev.off()
    }

  write.tree(phylo_tree, file = "snp_tree.nwk")
}

makePCAplot <- function(D) {
  png("./pca.png", width = 700, height = 700)
  pco1 <- dudi.pco(D, scannf=FALSE, nf=2)
  s.label(pco1$li*1.1, clab=0, pch="")
  colorplot(pco1$li, pco1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
  abline(v=0,h=0,col="grey", lty=2)
  title("Principal Coordinate Analysis\n-based on SNP distances-")
  dev.off()
  graphics.off()
}

main <- function(){
  ## Main function
  opt <- parseArgs()
  allele_matrix <- opt$allele_matrix

  D<-readAlleleMatrix(allele_matrix)
  mp2(D)

  makePCAplot(D)

}

main()
