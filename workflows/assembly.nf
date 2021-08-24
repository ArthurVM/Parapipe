// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {spades} from '../modules/assemblyModules.nf'
include {quast} from '../modules/assemblyModules.nf'
include {indexAssembly} from '../modules/assemblyModules.nf'
include {mapBT2} from '../modules/assemblyModules.nf'
include {pilon} from '../modules/assemblyModules.nf'

// define workflow
workflow assembly {

    take:
      input_files
      trimmed_fqs
      refdata

    main:
      spades(input_files, trimmed_fqs)

      quast(input_files, refdata, spades.out.scaffolds)

      indexAssembly(input_files, spades.out.scaffolds)

      mapBT2(input_files, trimmed_fqs, indexAssembly.out.bt2_index)

      pilon(input_files, mapBT2.out.bam, spades.out.scaffolds)

    emit:
      pilon_fasta = pilon.out.pilon_fasta
}
