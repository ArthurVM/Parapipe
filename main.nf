#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import modules
include {printHelp} from './modules/help.nf'
include {prepRef} from './workflows/prepref.nf'
include {preprocessing} from './workflows/preprocessing.nf'
include {callSNPs} from './workflows/varanalysis.nf'
include {callSTRs} from './workflows/varanalysis.nf'
include {assembly} from './workflows/assembly.nf'

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
if ( params.genome == "" ) {
    exit 1, "error: please provide a --genome argument. Use --help for an explenation for the parameters."
}
if ( params.pattern == "" ) {
    exit 1, "error: please provide a --pattern argument. Use --help for an explenation for the parameters."
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
    prepRef(params.genome)

    preprocessing(input_files, prepRef.out.ref_bt2index)

    callSNPs(input_files, preprocessing.out.bam, prepRef.out.refdata)

    assembly(input_files, preprocessing.out.trimmed_fqs, prepRef.out.refdata)

    callSTRs(input_files, preprocessing.out.trimmed_fqs, assembly.out.fasta, prepRef.out.refdata)
}
