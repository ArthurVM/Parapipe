// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../modules/preprocessingModules.nf'
include {countReads} from '../modules/preprocessingModules.nf'
include {fastp} from '../modules/preprocessingModules.nf'
include {fastQC} from '../modules/preprocessingModules.nf'
include {multiQC} from '../modules/preprocessingModules.nf'
include {trimGalore} from '../modules/preprocessingModules.nf'
include {map2Ref} from '../modules/preprocessingModules.nf'
include {picard} from '../modules/preprocessingModules.nf'
include {summarise} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
      ref_bt2index
      read_n_threshold

    main:
      // check fastq files for quality and read count
      checkFqValidity(input_files)
      countReads(checkFqValidity.out.checkValidity_fqs, read_n_threshold)
      fastp(countReads.out.countReads_fqs, read_n_threshold)

      // produce QC reports
      fastQC(fastp.out.fastp_fqs)
      fq_reports = fastQC.out.fastQC_report.collect()
      multiQC(fq_reports)

      // trim reads
      // temporarily removed
      // trimGalore(fastp.out.fastp_fqs)

      // map reads to reference sequence to produce BAM
      map2Ref(fastp.out.fastp_fqs, ref_bt2index)

      // run deduplication and read grouping on BAM
      picard(map2Ref.out.bam)

      // calculate mapping stats and capture in JSON
      summarise(picard.out.grouped_bam)

    emit:
      bam_pre = summarise.out.bam_pre_out
}
