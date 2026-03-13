// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {checkFqValidity} from '../modules/preprocessingModules.nf'
include {countReads} from '../modules/preprocessingModules.nf'
include {fastp} from '../modules/preprocessingModules.nf'
include {fastQC} from '../modules/preprocessingModules.nf'
include {multiQC} from '../modules/preprocessingModules.nf'
include {trimGalore} from '../modules/preprocessingModules.nf'

// define workflow
workflow preprocessing {

    take:
      input_files
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
      // DEPRECATED: handled by fastp
      // trimGalore(fastp.out.fastp_fqs)

    emit:
      cleaned_fqs = fastp.out.fastp_fqs
      multiQC_report = multiQC.out.multiQC_report
}
