#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import modules
include {printHelp} from './modules/help.nf'
include {prepRef} from './workflows/prepref.nf'
include {preprocessing} from './workflows/preprocessing.nf'
include {snp_analysis} from './workflows/snp_analysis.nf'

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

if ( params.ref == null ) {
    exit 1, "error: please provide a --ref argument. Use --help for an explenation for the parameters."
}

if ( params.pattern == null ) {
    exit 1, "error: please provide a --pattern argument. Use --help for an explenation for the parameters."
}

if ( params.yaml == null ) {
    exit 1, "error: please provide a --yaml argument. Use --help for an explenation for the parameters."
}
else {
  yaml = "${workflow.launchDir}/${params.yaml}"
}

if ( params.database == null ) {
  db = "false"
}
else {
  db = "${workflow.launchDir}/${params.database}"
}

log.info """
========================================================================
P A R A P I P E

Parameters used:
------------------------------------------------------------------------
--input_dir  ${params.input_dir}
--output_dir ${params.output_dir}
--ref		     ${params.ref}
--pattern		 ${params.pattern}
--yaml       ${params.yaml}
--database   ${db}

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
    /*******************************
    *    PREPREF WORKFLOW START    *
    ********************************/
    prepRef(params.ref)

    /*******************************
    * PREPROCESSING WORKFLOW START *
    ********************************/
    preprocessing(input_files, prepRef.out.ref_bt2index)

    /*******************************
    *      SNP WORKFLOW START      *
    ********************************/
    snp_analysis(db, preprocessing.out.bam, prepRef.out.refdata, preprocessing.out.mapstats_json, params.ref, yaml)

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
