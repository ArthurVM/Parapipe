// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {getMOI} from '../modules/phyloModules.nf'
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
      yaml

    main:
      // Investigate sample heterozygosity
      getMOI(bam_pre, refdata, ref_id)

      // Extract just bam files from tuple for collection
      formatPhyloInput(getMOI.out.bam_pp_moi)

      // Carry out SNP typing analysis
      wgSNV_phylo(formatPhyloInput.out.mapstats_json.collect(), formatPhyloInput.out.vcf.collect(), refdata, database)

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
      makeJSON(getMOI.out.bam_pp_moi, typingreport_json, wgSNV_phylo.out.al_stats_json)

      // Make report PDF for each sample
      makeSampleReports(getMOI.out.bam_pp_moi, env_json, makeJSON.out.report_json, wgSNV_phylo.out.snp_png)

      // Make report PDF for this run
      makeRunReport(env_json, makeJSON.out.report_json.collect(), wgSNV_phylo.out.snp_png, formatPhyloInput.out.moi_json.collect(), formatPhyloInput.out.moi_png.collect())

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
