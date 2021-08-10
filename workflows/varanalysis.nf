// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {runSNP} from '../modules/varanalysisModules.nf'

// define workflow
workflow callVariants {

    take:
      input_files
      bam
      refdata

    main:
      runSNP(input_files, bam, refdata)

}
