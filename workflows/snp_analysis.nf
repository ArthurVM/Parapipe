// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {callVariants} from '../modules/phyloModules.nf'
include {moi} from '../modules/phyloModules.nf'
include {molTyping_phylo} from '../modules/phyloModules.nf'
include {wgSNV_phylo} from '../modules/phyloModules.nf'
include {makeJSON} from '../modules/phyloModules.nf'
include {makeSampleReports} from '../modules/phyloModules.nf'
include {makeRunReport} from '../modules/phyloModules.nf'

// define SNP workflow
workflow phylo {

    take:
      input_files
      env_json
      database
      bam_pre
      refdata
      ref_id
      multiQC_report
      yaml
      mincov
    

    main:
      // Investigate sample heterogeneity
      callVariants(bam_pre, refdata, ref_id)

      // Carry out population Fws analysis
      moi(callVariants.out.vcf.collect(), refdata)

      // Extract just bam files from tuple for collection
      formatPhyloInput(callVariants.out.bam_vcf_tup)

      // Carry out SNP typing analysis
      wgSNV_phylo(formatPhyloInput.out.mapstats_json.collect(), formatPhyloInput.out.vcf.collect(), refdata, database, mincov)

      // if a yaml is provided then run in silico molecular typing
      if ( yaml != "false" ) {
        // Carry out molecular typing analysis
        molTyping_phylo(formatPhyloInput.out.bam.collect(), yaml)
        typingreport_json = molTyping_phylo.out.typingreport_json
      }

      // else capture that for the report JSON
      else {
        typingreport_json = "/None/"
      }

      // Make report JSON
      makeJSON(callVariants.out.bam_vcf_tup, typingreport_json, wgSNV_phylo.out.al_stats_json)

      // Make report PDF for each sample
      makeSampleReports(callVariants.out.bam_vcf_tup, env_json, makeJSON.out.report_json, wgSNV_phylo.out.snp_png, moi.out.moi_json, moi.out.moi_png)

      // Make report PDF for this run
      makeRunReport(env_json, makeJSON.out.report_json.collect(), wgSNV_phylo.out.snp_png, moi.out.moi_json, moi.out.moi_png, multiQC_report, wgSNV_phylo.out.dist_matrix)

}

process formatPhyloInput {

    input:
    tuple val(sample_name), path(bam), path(mapstats_jsons), path(vcf)

    output:
    path(bam), emit: bam
    path(mapstats_jsons), emit: mapstats_json
    path(vcf), emit: vcf

    script:
    """
    echo /${sample_name}/
    """
}
