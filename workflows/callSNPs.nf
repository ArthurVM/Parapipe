// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {runSNP} from '../modules/varanalysisModules.nf'
include {plotSNP} from '../modules/varanalysisModules.nf'
include {preprocessForPlotting} from '../modules/varanalysisModules.nf'
include {findSTRs} from '../modules/varanalysisModules.nf'

// define workflow
workflow callSNPs {

    take:
      input_files
      bam
      refdata
      pilon_fasta

    main:
      runSNP(input_files, bam, refdata)

      preprocessForPlotting(input_files, runSNP.out.vcf, refdata)

      plotSNP(input_files, preprocessForPlotting.out.vcf_path, preprocessForPlotting.out.preprocessed_fasta, refdata)

      findSTRs(input_files, pilon_fasta)
}