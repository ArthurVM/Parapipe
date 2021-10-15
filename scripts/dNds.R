library(ape)
library(optparse)

parseArgs <- function() {
  option_list = list(
    make_option(c("-a", "--aln"), default=NULL, help="alignment file in FASTA format.")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  if ( is.null(opt$aln)
  )
  {
    print_help(opt_parser)
    stop("All arguments need to be supplied.\n", call.=FALSE)
  }
  
  return(opt)
}

calcdNdS <- function(aln_file) {
  ## takes an alignment file in FASTA format and calculates dNdS for this region.
  aln <- read.FASTA(aln_file, type="DNA")
  print(aln)
  res <- dnds(aln, code=1, codonstart=1, quiet=FALSE)
  
  print(res)
}

main <- function() {
  # opt <- parseArgs()
  # calcdNdS(opt$aln)
  aln_file<-"/home/amorris/BioInf/Parapipe/test_wetrun/assdata/tmp.CHUDEA8_5000.aln.fa"
  calcdNdS(aln_file)
}

main()