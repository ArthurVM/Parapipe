// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {findSNPs} from '../../modules/snpAnalysisModules.nf'
include {plotSNP} from '../../modules/snpAnalysisModules.nf'
include {preprocessForPlotting} from '../../modules/snpAnalysisModules.nf'
include {makeVCFtree} from '../../modules/snpAnalysisModules.nf'
include {makeJSON} from '../../modules/snpAnalysisModules.nf'
include {makeSNPreport} from '../../modules/snpAnalysisModules.nf'

// define SNP workflow
workflow SNP_analysis {

    take:
      input_files
      database
      bam
      refdata
      mapstats_json
      ref_id

    main:
      findSNPs(input_files, bam, refdata, ref_id)

      makeVCFtree(findSNPs.out.gvcf.collect(), database, refdata, findSNPs.out.merge_and_plot)

      preprocessForPlotting(input_files, findSNPs.out.vcf, refdata)

      plotSNP(input_files, preprocessForPlotting.out.vcf_path, preprocessForPlotting.out.fasta_path, refdata)

      makeJSON(input_files, mapstats_json, makeVCFtree.out.newick_tree)

      makeSNPreport(input_files, makeJSON.out.snpreport_json, plotSNP.out.snpQC_plots)

}
