// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {getMOI} from '../modules/snpAnalysisModules.nf'
include {phylo} from '../modules/snpAnalysisModules.nf'
include {makeJSON} from '../modules/snpAnalysisModules.nf'
include {makeSampleReports} from '../modules/snpAnalysisModules.nf'
include {makeRunReport} from '../modules/snpAnalysisModules.nf'

// define SNP workflow
workflow snp_analysis {

    take:
      input_files
      env_json
      database
      bam_pre
      refdata
      ref_id
      yaml

    main:
      // Investigate sample heterozygosity
      getMOI(bam_pre, refdata, ref_id)

      // Extract just bam files from tuple for collection
      formatPhyloInput(getMOI.out.bam_pp_moi)

      // Carry out SNP and ST phylogenetic analysis
      phylo(formatPhyloInput.out.bam.collect(), formatPhyloInput.out.mapstats_json.collect(), formatPhyloInput.out.vcf.collect(), refdata, database, yaml)

      // Make report JSON
      makeJSON(getMOI.out.bam_pp_moi, phylo.out.snp_nwk, phylo.out.iqtree_json, phylo.out.st_stats, phylo.out.al_stats_json)

      // Make report PDF for each sample
      makeSampleReports(getMOI.out.bam_pp_moi, env_json, makeJSON.out.report_json, phylo.out.snp_png)

      // Make report PDF for this run
      makeRunReport(env_json, makeJSON.out.report_json.collect(), phylo.out.snp_png, formatPhyloInput.out.moi_json.collect(), formatPhyloInput.out.moi_png.collect())

}

process formatPhyloInput {

    input:
    tuple val(sample_name), path(bam), path(mapstats_jsons), path(moi_json), path(moi_png), path(vcf)

    output:
    path(bam), emit: bam
    path(mapstats_jsons), emit: mapstats_json
    path(moi_json), emit: moi_json
    path(moi_png), emit: moi_png
    path(vcf), emit: vcf

    script:
    """
    echo /${sample_name}/
    """
}
