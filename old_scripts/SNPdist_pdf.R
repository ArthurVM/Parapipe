#!/usr/bin/env Rscript

# version 1, 19 AUG 21, A.V. Morris

library(vcfR)
library(ape)
library(optparse)
library(glue)

parseArgs <- function() {
  ## Simple argument parser
  
    option_list = list(
        make_option(c("-v", "--vcf_dir"), type="character", default=NULL, 
              help="directory containing vcf files split by chromosome", metavar="character"),
        make_option(c("-f", "--fasta_file"), type="character", default=NULL,
              help="fasta file", metavar="character"),
        make_option(c("-g", "--gff_file"), type="character", default=NULL,
              help="gff3 file path", metavar="character")
    ); 

    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if ( is.null(opt$vcf_dir) | is.null(opt$fasta_file) | is.null(opt$gff_file) ){
         print_help(opt_parser)
         stop("All arguements need to be supplied.\n", call.=FALSE)
    }
    
    return(opt)
}

readFiles <- function(opt){
  ## Reads files into memory and subsets by chromosome ID
  
  fasta_dna <- ape::read.dna(opt$fasta_file, format = "fasta", comment.char = "|")
  gff       <- read.table(opt$gff_file, sep="\t", quote="")
  
  pdf("SNPDist.pdf", width = 10, height = 10)
  for (seqid in names(fasta_dna)){
    vcf_file   <- glue("{opt$vcf_dir}/tmp.{seqid}.vcf")
    vcf       <- read.vcfR(vcf_file, verbose = TRUE )
    
    fastasubset <- fasta_dna[ grep( seqid, names(fasta_dna) ) ]
    gffsubset <- gff[gff$V1 %in% seqid, ]
    
    plotter(fastasubset, gffsubset, vcf, seqid)
  }
  dev.off()
}
plotter <- function(fastasubset, gffsubset, vcf, seqid){
  ## generate plots
  
  chrome <- create.chromR(name=paste(seqid), vcf=vcf, seq=fastasubset, ann=gffsubset)
  chrome.nomasker <- proc.chromR(chrome, verbose=TRUE)
  
  plot(chrome.nomasker,main=seqid)
  chromoqc(chrome.nomasker, dp.alpha=40)
}

main <- function(){
  ## Main function
  opt <- parseArgs()
  readFiles(opt)
}

main()

