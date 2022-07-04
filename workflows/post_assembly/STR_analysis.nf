// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {findSTRs} from '../../modules/varanalysisModules.nf'
include {indexScaffolds} from '../../modules/varanalysisModules.nf'
include {map2Scaffolds} from '../../modules/varanalysisModules.nf'
include {runSTRViper} from '../../modules/varanalysisModules.nf'
include {runSNPEff} from '../../modules/varanalysisModules.nf'

// define STR workflow
workflow STR_analysis {

    take:
      input_files
      trimmed_fqs
      scaffolds
      refdata

    main:
      findSTRs(input_files, scaffolds)

      indexScaffolds(input_files, scaffolds)

      map2Scaffolds(input_files, trimmed_fqs, indexScaffolds.out.scaffold_bt2index)

      runSTRViper(input_files, map2Scaffolds.out.bam, findSTRs.out.trf_dat, scaffolds)

      // SNPEff(input_files, STRViper.out.str_vcf, refdata)
}
