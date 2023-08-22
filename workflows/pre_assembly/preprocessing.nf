// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../../modules/preprocessingModules.nf'
include {countReads} from '../../modules/preprocessingModules.nf'
include {fastp} from '../../modules/preprocessingModules.nf'
include {fastQC} from '../../modules/preprocessingModules.nf'
include {multiQC} from '../../modules/preprocessingModules.nf'
include {trimGalore} from '../../modules/preprocessingModules.nf'
include {map2Ref} from '../../modules/preprocessingModules.nf'
include {picard} from '../../modules/preprocessingModules.nf'
include {gini} from '../../modules/preprocessingModules.nf'
include {summarise} from '../../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
      ref_bt2index

    main:
      checkFqValidity(input_files)

      countReads(checkFqValidity.out.checkValidity_fqs)

      fastp(countReads.out.countReads_fqs)

      fastQC(fastp.out.fastp_fqs)

      fq_reports = fastQC.out.fastQC_report.collect()

      multiQC(fq_reports)

      // trimGalore(fastp.out.fastp_fqs)

      map2Ref(fastp.out.fastp_fqs, ref_bt2index)

      picard(checkFqValidity.out.checkValidity_fqs, map2Ref.out.bam)

      gini(fastp.out.fastp_fqs, picard.out.grouped_bam)

      summarise(fastp.out.fastp_fqs, picard.out.grouped_bam)

    emit:
      bam = picard.out.grouped_bam
      trimmed_fqs = fastp.out.fastp_fqs
      mapstats_json = summarise.out.mapstats_json
}
