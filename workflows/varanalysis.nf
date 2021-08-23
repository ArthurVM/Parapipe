// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {runSNP} from '../modules/varanalysisModules.nf'
include {plotSNP} from '../modules/varanalysisModules.nf'
include {preprocessForPlotting} from '../modules/varanalysisModules.nf'

// define workflow
workflow callVariants {

    take:
      input_files
      bam
      refdata

    main:
      runSNP(input_files, bam, refdata)

      preprocessForPlotting(input_files, runSNP.out.vcf, refdata)

      plotSNP(input_files, preprocessForPlotting.out.vcf_path, preprocessForPlotting.out.preprocessed_fasta, refdata)
}
