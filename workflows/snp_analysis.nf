// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {getMOI} from '../modules/snpAnalysisModules.nf'
include {phylo} from '../modules/snpAnalysisModules.nf'
// include {findSNPs} from '../modules/snpAnalysisModules.nf'
// include {makeVCFtree} from '../modules/snpAnalysisModules.nf'
include {makeJSON} from '../modules/snpAnalysisModules.nf'
include {makeReport} from '../modules/snpAnalysisModules.nf'

// define SNP workflow
workflow snp_analysis {

    take:
      database
      bam
      refdata
      mapstats_json
      ref_id
      yaml

    main:
      getMOI(bam, refdata, ref_id)

      formatPhyloInput(bam)

      phylo(formatPhyloInput.out.bam.collect(), refdata, database, yaml)

      // findSNPs(bam, refdata, ref_id)
      //
      // makeVCFtree(findSNPs.out.gvcf.collect(), database, refdata, findSNPs.out.merge_and_plot)

      makeJSON(bam, mapstats_json, phylo.out.snp_nwk, phylo.out.iqtree_json)

      // makeSNPreport(input_files, makeJSON.out.snpreport_json, plotSNP.out.snpQC_plots)

}

process formatPhyloInput {

    input:
    tuple val(sample_name), path(bam)

    output:
    path(bam), emit: bam

    script:
    """
    echo /${sample_name}/
    """
}
