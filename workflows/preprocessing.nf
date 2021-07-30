// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {fastQC} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files

    main:
      fastQC(input_files)
}
