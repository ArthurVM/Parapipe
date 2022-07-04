#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import modules
include {printHelp} from './modules/help.nf'
include {prepRef} from './workflows/pre_assembly/prepref.nf'
include {preprocessing} from './workflows/pre_assembly/preprocessing.nf'
include {SNP_analysis} from './workflows/pre_assembly/analyse_variants.nf'
include {targetAnalysis} from './workflows/pre_assembly/analyse_variants.nf'
include {assembly} from './workflows/post_assembly/assembly.nf'
include {STR_analysis} from './workflows/post_assembly/STR_analysis.nf'
include {postAssemblyAnalysis} from './workflows/post_assembly/postAssemblyAnalysis.nf'

/*
 ANSI escape codes to allow colour-coded output messages
 This code is from https://github.com/angelovangel
 */

ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

params.help = ""

if (params.help) {
    printHelp()
    exit(0)
}

// check mandatory parameters
if ( params.input_dir == null ) {
    exit 1, "error: please provide an --input_dir argument. Use --help for an explenation for the parameters."
}

if ( params.output_dir == null ) {
    exit 1, "error: please provide an --output_dir argument. Use --help for an explenation for the parameters."
}

if ( params.genome == null ) {
    exit 1, "error: please provide a --genome argument. Use --help for an explenation for the parameters."
}

if ( params.pattern == null ) {
    exit 1, "error: please provide a --pattern argument. Use --help for an explenation for the parameters."
}

if ( params.only_SNP == null ) {
    params.only_SNP = false
}

// if ( ( params.only_SNP != "yes" ) && ( params.only_SNP != "no" ) ) {
//     exit 1, "error: --only_SNP is mandatory and must be either \"yes\" or \"no\"."
// }

if ( params.ref_scaffold == null ) {
    params.ref_scaffold = false
}

// if ( ( params.ref_scaffold != "yes" ) && ( params.ref_scaffold != "no" ) ) {
//     exit 1, "error: --ref_scaffold is mandatory and must be either \"yes\" or \"no\"."
// }

if ( ( params.only_SNP == true ) && ( params.ref_scaffold == true ) ) {
    exit 1, "error: --ref_scaffold and --only_SNP are mutually exclusive."
}

if ( params.target_flanks == null ) {
    params.target_flanks = false
}

log.info """
========================================================================
P A R A P I P E

Parameters used:
------------------------------------------------------------------------
--input_dir ${params.input_dir}
--output_dir ${params.output_dir}
--genome		${params.genome}
--pattern		${params.pattern}
--ref_scaffold  ${params.ref_scaffold}
--only_SNP  ${params.only_SNP}
--target_flanks ${params.target_flanks}

Runtime data:
------------------------------------------------------------------------
Profile     ${workflow.profile}
User        ${workflow.userName}
Launch dir  ${workflow.launchDir}
"""
.stripIndent()

// main pipeline workflow
workflow {

  pattern = params.pattern
	reads = params.input_dir + pattern
	numfiles = file(reads) // count the number of files

  Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
   .ifEmpty { error "cannot find any reads matching ${pattern} in ${indir}" }
   .set{ input_files }

  main:
    ref_scaffold_bool = params.ref_scaffold
    only_SNP_bool = params.only_SNP
    prepRef(params.genome)

    preprocessing(input_files, prepRef.out.ref_bt2index)

    SNP_analysis(input_files, preprocessing.out.bam, prepRef.out.refdata, preprocessing.out.mapstats_json)

    if ( params.target_flanks != false ) {
      targetAnalysis(input_files, preprocessing.out.trimmed_fqs, params.target_flanks)
    }

    // check whether to run the rest of the pipeline or stop after calling SNPs
    if ( only_SNP_bool == false ) {
      // run the whole pipeline
      assembly(input_files, preprocessing.out.trimmed_fqs, ref_scaffold_bool, prepRef.out.refdata)

      assemblies_list = assembly.out.fasta.collect()
      annotations_list = assembly.out.gff.collect()

      /* TODO: work on this block
      * if ( ref_scaffold_bool == "yes" ) {
      *   postAssemblyAnalysis(assemblies_list, annotations_list, prepRef.out.refdata)
      * }
      */

      STR_analysis(input_files, preprocessing.out.trimmed_fqs, assembly.out.fasta, prepRef.out.refdata)
    }
}

workflow.onComplete {
    if ( workflow.success ) {
        log.info """
        ============================================
        ${ANSI_GREEN}PARAPIPE completed successfully
        """
        .stripIndent()
    }
    else {
        log.info """
        ===========================================
        ${ANSI_RED}Finished with errors${ANSI_RESET}
        """
        .stripIndent()
    }
}
