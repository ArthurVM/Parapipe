// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {findSNPs} from '../modules/snpAnalysisModules.nf'
include {getMOI} from '../modules/snpAnalysisModules.nf'
include {makeVCFtree} from '../modules/snpAnalysisModules.nf'
include {makeJSON} from '../modules/snpAnalysisModules.nf'
include {makeSNPreport} from '../modules/snpAnalysisModules.nf'

// define SNP workflow
workflow snp_analysis {

    take:
      database
      bam
      refdata
      mapstats_json
      ref_id

    main:
      getMOI(bam, refdata, ref_id)

      findSNPs(bam, refdata, ref_id)

      makeVCFtree(findSNPs.out.gvcf.collect(), database, refdata, findSNPs.out.merge_and_plot)

      makeJSON(bam, mapstats_json, makeVCFtree.out.newick_tree)

      // makeSNPreport(input_files, makeJSON.out.snpreport_json, plotSNP.out.snpQC_plots)

}
