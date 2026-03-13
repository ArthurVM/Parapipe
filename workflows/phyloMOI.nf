// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {snippyCall} from '../modules/snippyModules.nf'
include {summarise} from '../modules/preprocessingModules.nf'
include {moi} from '../modules/phyloModules.nf'
include {molTyping_phylo} from '../modules/phyloModules.nf'
include {wgSNV_phylo} from '../modules/phyloModules.nf'
include {gp60BlooMine} from '../modules/phyloModules.nf'
include {makeJSON} from '../modules/phyloModules.nf'
include {makeSampleReports} from '../modules/phyloModules.nf'
include {makeRunReport} from '../modules/phyloModules.nf'

// define SNP workflow
workflow phyloMOI {

    take:
      env_json
      database
      cleaned_fqs
      refdata
      gp60_bt2_db
      ref_id
      multiQC_report
      yaml
      mincov
      missing
      maf
      mac
      pattern
    

    main:
      // Map reads and call variants with Snippy
      // fastp emits stdout as a 4th tuple element; drop it for Snippy
      snippy_inputs = cleaned_fqs.map{ sample, fq1, fq2, _ -> tuple(sample, fq1, fq2) }
      snippyCall(snippy_inputs, refdata, ref_id)

      // Generate mapping stats from Snippy BAMs
      summarise(snippyCall.out.bam_vcf.map{ sample_name, bam, vcf -> tuple(sample_name, bam) })

      // Carry out population Fws analysis
      moi(snippyCall.out.raw_vcf.collect(), refdata, missing, maf, mac)

      // Build combined channel for downstream reporting
      snippy_reports = summarise.out.bam_pre_out.join(snippyCall.out.bam_vcf.map{ sample_name, bam, vcf -> tuple(sample_name, vcf) })

      // gp60 subtype analysis (Cryptosporidium only)
      def gp60_reports
      if (ref_id.toLowerCase().startsWith('cryptosporidium')) {
        def gp60Probes = file("${baseDir}/resources/gp60_probes/gp60_probes.fa")
        gp60BlooMine(snippy_inputs, pattern, gp60Probes, gp60_bt2_db)
        gp60_reports = gp60BlooMine.out.filtered_gp60_report
      }
      else {
        gp60_reports = snippy_reports.map{ sample_name, bam, mapstats, vcf -> tuple(sample_name, "None") }
      }

      snippy_reports_gp60 = snippy_reports.join(gp60_reports)

      // Carry out SNP typing analysis
      wgSNV_phylo(summarise.out.bam_pre_out.map{ sample, bam, mapstats -> mapstats }.collect(), snippyCall.out.vcf.collect(), refdata, database, mincov)

      // // Make report JSON
      makeJSON(snippy_reports_gp60, wgSNV_phylo.out.al_stats_json)

      // // Make report for each sample
      // makeSampleReports(snippy_reports, env_json, makeJSON.out.report_json, wgSNV_phylo.out.snp_png, moi.out.moi_json, moi.out.moi_png)

      // // Make report for this run
      makeRunReport(env_json, makeJSON.out.report_json.collect(), wgSNV_phylo.out.snp_png, moi.out.moi_json, moi.out.moi_png, multiQC_report, wgSNV_phylo.out.dist_matrix)

}
