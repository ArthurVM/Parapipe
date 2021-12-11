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
include {summarise} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
      ref_bt2index

    main:
      checkFqValidity(input_files)

      fastQC(checkFqValidity.out.checkValidity_fqs)

      fq_reports = fastQC.out.fastQC_report.collect()

      multiQC(fq_reports)

      trimGalore(checkFqValidity.out.checkValidity_fqs)

      map2Ref(checkFqValidity.out.checkValidity_fqs, trimGalore.out.tg_fqs, ref_bt2index)

      // deduplication(checkFqValidity.out.checkValidity_fqs, map2Ref.out.bam)

      gini(checkFqValidity.out.checkValidity_fqs, map2Ref.out.bam)

      summarise(checkFqValidity.out.checkValidity_fqs, map2Ref.out.bam, gini.out.GG)

    emit:
      bam = map2Ref.out.bam
      trimmed_fqs = trimGalore.out.tg_fqs
}
