// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {fastQC} from '../modules/preprocessingModules.nf'
include {multiQC} from '../modules/preprocessingModules.nf'
include {checkFqValidity} from '../modules/preprocessingModules.nf'
include {trimGalore} from '../modules/preprocessingModules.nf'
include {map2Ref} from '../modules/preprocessingModules.nf'
include {deduplication} from '../modules/preprocessingModules.nf'
include {gini} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
      ref_bt2index

    main:
      checkFqValidity(input_files)

      fastQC(input_files)

      // multiQC(fastQC.out.fastQC_report.toList())

      trimGalore(input_files)

      map2Ref(input_files, trimGalore.out.tg_fqs, ref_bt2index)

      deduplication(input_files, map2Ref.out.bam)

      gini(input_files, map2Ref.out.bam)

    emit:
      bam = map2Ref.out.bam
      trimmed_fqs = trimGalore.out.tg_fqs
}
