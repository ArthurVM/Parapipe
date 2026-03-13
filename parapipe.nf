#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import modules
include {printHelp} from './modules/help.nf'
include {captureEnv} from './modules/captureEnv.nf'
include {prepRef} from './workflows/prepref.nf'
include {preprocessing} from './workflows/preprocessing.nf'
include {phyloMOI} from './workflows/phyloMOI.nf'

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
if ( params.input_dir == "" ) {
    exit 1, "error: please provide an --input_dir argument. Use --help for an explenation for the parameters."
}

if ( params.output_dir == "" ) {
    exit 1, "error: please provide an --output_dir argument. Use --help for an explenation for the parameters."
}

if ( params.ref == "" ) {
    exit 1, "error: please provide a --ref argument. Use --help for an explenation for the parameters."
}

if ( params.pattern == "" ) {
    exit 1, "error: please provide a --pattern argument. Use --help for an explenation for the parameters."
}

// initialise parameters
if ( params.yaml == "false" ) {
    yaml = params.yaml
}
else {
    print("WARNING: typing functionality will be deprecated in version 0.2")
    yaml = params.yaml
}

read_n_threshold = params.read_n_threshold
database = params.database
mincov = params.mincov
missing = params.missing
maf = params.maf
mac = params.mac

log.info """
===================================================================================================
                  ██████╗░░█████╗░██████╗░░█████╗░██████╗░██╗██████╗░███████╗
                  ██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔══██╗██║██╔══██╗██╔════╝
                  ██████╔╝███████║██████╔╝███████║██████╔╝██║██████╔╝█████╗░░
                  ██╔═══╝░██╔══██║██╔══██╗██╔══██║██╔═══╝░██║██╔═══╝░██╔══╝░░
                  ██║░░░░░██║░░██║██║░░██║██║░░██║██║░░░░░██║██║░░░░░███████╗
                  ╚═╝░░░░░╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░╚═╝╚═╝░░░░░╚═╝╚═╝░░░░░╚══════╝
===================================================================================================

Parameters used:
---------------------------------------------------------------------------------------------------
--input_dir             ${params.input_dir}
--output_dir            ${params.output_dir}
--ref                   ${params.ref}
--pattern               ${params.pattern}
--read_n_threshold      ${read_n_threshold}
--yaml                  ${yaml}
--database              ${database}
--mincov                ${mincov}
--missing               ${missing}
--maf                   ${maf}
--mac                   ${mac}

Runtime data:
---------------------------------------------------------------------------------------------------
Profile                 ${workflow.profile}
User                    ${workflow.userName}
Launch dir              ${workflow.launchDir}
"""
.stripIndent()

// main pipeline workflow
workflow {

    pattern = params.pattern
	reads = params.input_dir + pattern
	numfiles = file(reads) // count the number of files

    Channel.fromFilePairs(reads, flat: true, checkIfExists: true, size: -1)
        .ifEmpty { error "cannot find any files matching ${pattern} in ${params.input_dir}" }
        .set{ input_files }

  main:

    /************************************************
    *               ENV CAPTURE START               *
    *************************************************/
    captureEnv( params.input_dir, \
                params.output_dir, \
                params.ref, \
                yaml, \
                database, \
                read_n_threshold, \
                mincov, \
                missing, \
                maf, \
                mac )



    /************************************************
    *            PREPREF WORKFLOW START             *
    *************************************************/
    prepRef(params.ref)



    /************************************************
    *          PREPROCESSING WORKFLOW START         *
    *************************************************/
    preprocessing(  input_files, \
                    read_n_threshold    )



    /************************************************
    *      PARAPIPE SNP AND MOI WORKFLOW START      *
    *************************************************/
    phyloMOI(   captureEnv.out.env_json, \
                database, \
                preprocessing.out.cleaned_fqs, \
                prepRef.out.refdata, \
                prepRef.out.gp60_bt2_db, \
                params.ref, \
                preprocessing.out.multiQC_report, \
                yaml, \
                mincov, \
                missing, \
                maf, \
                mac, \
                pattern )
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
