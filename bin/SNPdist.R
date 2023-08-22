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
        make_option(c("-f", "--fasta_dir"), type="character", default=NULL,
              help="fasta dir", metavar="character"),
        make_option(c("-g", "--gff_file"), type="character", default=NULL,
              help="gff3 file path", metavar="character")
    );

    opt_parser = OptionParser(option_list=option_list);
    opt = parse_args(opt_parser);

    if ( is.null(opt$vcf_dir) | is.null(opt$fasta_dir) | is.null(opt$gff_file) ){
         print_help(opt_parser)
         stop("All arguements need to be supplied.\n", call.=FALSE)
    }

    return(opt)
}

readFiles <- function(fasta_dir, gff_file, vcf_dir){
  ## Reads files into memory and subsets by chromosome ID

  for (fasta_file in list.files(fasta_dir)) {
    prefix <- strsplit(fasta_file, ".fasta")[[1]]
    seqid <- strsplit(prefix, ".tmp")[[1]]

    fasta_path <- glue("{fasta_dir}/{fasta_file}")
    vcf_file  <- glue("{vcf_dir}/{prefix}.vcf")

    print(fasta_path, vcf_file)

    if ( file.exists(fasta_path) && file.exists(vcf_file) ) {
      fasta_dna  <- ape::read.dna(fasta_path, format = "fasta", comment.char = "|")

      vcf       <- read.vcfR(vcf_file, verbose = TRUE )

      gff       <- read.table(gff_file, sep="\t", quote="")
      gffsubset <- gff[gff$V1 %in% seqid, ]
     
      plotter(fasta_dna, gffsubset, vcf, seqid)
    }
  }
}

plotter <- function(fastasubset, gffsubset, vcf, seqid){
  ## generate plots

  chrome <- create.chromR(name=paste(seqid), vcf=vcf, seq=fastasubset, ann=gffsubset)
  chrome.nomasker <- proc.chromR(chrome, verbose=TRUE)

  jpeg(glue("{seqid}_chromQC.jpeg"), width = 1000, height = 1000)
  chromoqc(chrome.nomasker, dp.alpha=40)
  dev.off()

  jpeg(glue("{seqid}_chromPlot.jpeg"), width = 1000, height = 1000)
  plot(chrome.nomasker,main=seqid)
  dev.off()
}

main <- function(){
  ## Main function
  # fasta_dir <- "/home/amorris/BioInf/Parapipe/work/76/c7b6b3385a8fb7f0f7808d5492d435/fasta_split/"
  # gff_file <- "/home/amorris/BioInf/Parapipe/work/76/c7b6b3385a8fb7f0f7808d5492d435/mycobacterium_interjectum.gff"
  # vcf_dir <- "/home/amorris/BioInf/Parapipe/work/76/c7b6b3385a8fb7f0f7808d5492d435/vcf_split/"
  # readFiles(fasta_dir, gff_file, vcf_dir)
  
  opt <- parseArgs()
  readFiles(opt$fasta_dir, opt$gff_file, opt$vcf_dir)
}

main()
