library("vcfR")
library("adegenet")
library("ape")
library("optparse")
library("glue")

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

makePhylo <- function(vcf) {
  vcf <- "/home/amorris/BioInf/Parapipe/test_wetrunBIG_parv/GVCFs/merged.vcf"
  gind <- vcfR2genind(read.vcfR( vcf , verbose = FALSE ))
  row.names(gind@tab)
  D <- dist(tab(gind))
  tre <- nj(D)
  for ( s in tre$"tip.label" ) {
    fout <- glue("/home/amorris/BioInf/Parapipe/test_wetrunBIG_parv/{s}.pdf")
    pdf(fout, width = 10, height = 10)
    plot(tre, type="unrooted", edge.w=2, tip.color = c("red")[(tre$"tip.label"!=s)+1])
    edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
    dev.off()
  }
}

main <- function(){
  ## Main function
  opt <- parseArgs()
  makePhylo(opt$vcf)
}

main()