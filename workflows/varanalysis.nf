// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {findSNPs} from '../modules/varanalysisModules.nf'
include {plotSNP} from '../modules/varanalysisModules.nf'
include {preprocessForPlotting} from '../modules/varanalysisModules.nf'
include {findSTRs} from '../modules/varanalysisModules.nf'
include {indexPilonFasta} from '../modules/varanalysisModules.nf'
include {map2PilonFasta} from '../modules/varanalysisModules.nf'
include {processSTRs} from '../modules/varanalysisModules.nf'

// define SNP workflow
workflow callSNPs {

    take:
      input_files
      bam
      refdata

    main:
      findSNPs(input_files, bam, refdata)

      preprocessForPlotting(input_files, findSNPs.out.vcf, refdata)

      plotSNP(input_files, preprocessForPlotting.out.vcf_path, preprocessForPlotting.out.preprocessed_fasta, refdata)
}

// define STR workflow
workflow callSTRs {

    take:
      input_files
      trimmed_fqs
      pilon_fasta

    main:
      findSTRs(input_files, pilon_fasta)

      indexPilonFasta(input_files, pilon_fasta)

      map2PilonFasta(input_files, trimmed_fqs, indexPilonFasta.out.bt2_pilon_index)

      processSTRs(input_files, map2PilonFasta.out.bam, findSTRs.out.trf_dat, pilon_fasta)
}
