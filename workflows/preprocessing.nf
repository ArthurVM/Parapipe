// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {fastQC} from '../modules/preprocessingModules.nf'
include {checkFqValidity} from '../modules/preprocessingModules.nf'
include {trimGalore} from '../modules/preprocessingModules.nf'
include {mapBT2} from '../modules/preprocessingModules.nf'
include {deduplication} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
      ref_bt2index

    main:
      checkFqValidity(input_files)
      fastQC(input_files)
      trimGalore(input_files)
      mapBT2(input_files, trimGalore.out.tg_fqs, ref_bt2index)
      deduplication(input_files, mapBT2.out.bam)
}
