// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {findSNPs} from '../../modules/snpAnalysisModules.nf'
include {plotSNP} from '../../modules/snpAnalysisModules.nf'
include {preprocessForPlotting} from '../../modules/snpAnalysisModules.nf'
include {makeVCFtree} from '../../modules/snpAnalysisModules.nf'
include {makeSNPreport} from '../../modules/snpAnalysisModules.nf'
include {runBlooMine} from '../../modules/targetLocusAnalysisModules.nf'

// define SNP workflow
workflow SNP_analysis {

    take:
      input_files
      bam
      refdata
      mapstats_json

    main:
      findSNPs(input_files, bam, refdata)

      makeVCFtree(findSNPs.out.gvcf.collect(), refdata, findSNPs.out.merge_and_plot)

      preprocessForPlotting(input_files, findSNPs.out.vcf, refdata)

      plotSNP(input_files, preprocessForPlotting.out.vcf_path, preprocessForPlotting.out.preprocessed_fasta, refdata)

      makeSNPreport(input_files, mapstats_json, makeVCFtree.out.newick_tree)

}

// define target analysis workflow
workflow targetAnalysis {

    take:
      input_files
      trimmed_fqs
      target_flanks

    main:
      runBlooMine(input_files, trimmed_fqs, target_flanks)
}
