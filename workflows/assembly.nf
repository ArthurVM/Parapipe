// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {spades} from '../modules/assemblyModules.nf'
include {quast} from '../modules/assemblyModules.nf'

// define workflow
workflow assembly {

    take:
      input_files
      trimmed_fqs
      refdata

    main:
      spades(input_files, trimmed_fqs)

      quast(input_files, refdata, spades.out.scaffolds)
}
