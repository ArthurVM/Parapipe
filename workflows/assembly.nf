// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {spades} from '../modules/assemblyModules.nf'
include {quast} from '../modules/assemblyModules.nf'
include {indexAssembly} from '../modules/assemblyModules.nf'
include {map2SPAdesFasta} from '../modules/assemblyModules.nf'
include {pilon} from '../modules/assemblyModules.nf'
include {abacas} from '../modules/assemblyModules.nf'
include {liftover} from '../modules/assemblyModules.nf'

// define workflow
workflow assembly {

    take:
      input_files
      trimmed_fqs
      ref_scaffold_bool
      refdata

    main:
      spades(input_files, trimmed_fqs)

      quast(input_files, refdata, spades.out.scaffolds)

      indexAssembly(input_files, spades.out.scaffolds)

      map2SPAdesFasta(input_files, trimmed_fqs, indexAssembly.out.bt2_index)

      pilon(input_files, map2SPAdesFasta.out.bam, spades.out.scaffolds)

      if ( ref_scaffold_bool == "yes" ) {
        // run ABACAS if reference guided scaffolding is requested
        abacas(input_files, pilon.out.pilon_fasta, refdata)
        scaffolds_fasta = abacas.out.abacas_fasta
      }
      else {
        scaffolds_fasta = pilon.out.pilon_fasta
      }

      liftover(input_files, scaffolds_fasta, refdata)

    emit:
      // pilon_fasta = pilon.out.pilon_fasta
      fasta = scaffolds_fasta
}
