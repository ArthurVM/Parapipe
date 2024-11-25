# version 1, 19 AUG 21, A.V. Morris

library(optparse)
library(glue)
library(SeqArray)


parseArgs <- function() {
  ## Simple argument parser
  
  option_list = list(
    make_option(c("-v", "--vcf"), type="character", default=NULL,
                help="path to the merged vcf file.", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  if ( is.null(opt$vcf) ){
    print_help(opt_parser)
    stop("All arguments need to be supplied.\n", call.=FALSE)
  }
  
  return(opt)
}


get_MAF <- function(gdsfile) {
  stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
  vars <- seqSummary(gdsfile, check="none", verbose=FALSE)$format$ID
  if(!("AD" %in% vars)) {
    stop("Must have annotaion/format/AD tag to compute minor allele frequencies")
  }
  
  maf <- function(x) {
    # frequency of ref and alt reads over all samples
    coverage <- colSums(x, na.rm = TRUE)
    min(coverage / sum(coverage))
  }
  res <- seqApply(gdsfile, "annotation/format/AD",
                  function(x) maf(x),
                  margin = "by.variant",
                  as.is = "double")
  return(res)
}


get_heterozygosity <- function (gdsfile) 
{
  maf <- get_MAF(gdsfile)
  return(1 - (maf^2 + (1 - maf)^2))
}


get_heterozygosity_by_sample <- function (gdsfile) 
{
  stopifnot(inherits(gdsfile, "SeqVarGDSClass"))
  vars <- seqSummary(gdsfile, check = "none", verbose = FALSE)$format$ID
  if (!("AD" %in% vars)) {
    stop("Must have annotaion/format/AD tag to compute B-allele frequencies")
  }
  heterozygosity <- function(x) {
    depth <- rowSums(x)
    p <- x[, 1]/depth
    q <- x[, 2]/depth
    1 - (p^2 + q^2)
  }
  het <- seqApply(gdsfile, "annotation/format/AD", heterozygosity, 
                  margin = "by.variant", as.is = "list")
  het <- matrix(unlist(het), ncol = length(het), dimnames = list(sample = seqGetData(gdsfile, 
                                                                                     "sample.id"), variant = seqGetData(gdsfile, "variant.id")))
  return(het)
}


get_Fws <- function (gdsfile) 
{
  sample.het <- get_heterozygosity_by_sample(gdsfile)
  population.het <- get_heterozygosity(gdsfile)
  maf.bins <- findInterval(get_MAF(gdsfile), seq(0, 0.5, length.out = 11))
  mu.population.het <- tapply(population.het, maf.bins, mean)
  mu.sample.het <- apply(sample.het, 1, function(x) tapply(x, 
                                                           maf.bins, mean, na.rm = TRUE))
  fws <- apply(mu.sample.het, 2, function(x) 1 - lm(x ~ mu.population.het - 
                                                      1)$coeff)
  return(fws)
}


main <- function() {
  
  opt <- parseArgs()
  prefix <- strsplit(opt$vcf, ".vcf")[[1]]
  gds_file <- glue("tmp_vcf.gds")
  print(gds_file)
  seqVCF2GDS(opt$vcf, gds_file)
  
  gds_handle <- seqOpen(gds_file)
  fws <- 1-get_Fws(gds_handle)
  fws_df <- as.data.frame(fws)
  
  write.csv(fws_df, "./fws.csv")
}

main()