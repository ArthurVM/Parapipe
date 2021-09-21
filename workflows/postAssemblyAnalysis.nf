// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {makeChromosomeFastas} from '../modules/postAssemblyAnalysisModules.nf'
include {chromosomeAlignment} from '../modules/postAssemblyAnalysisModules.nf'

// define workflow
workflow postAssemblyAnalysis {

    take:
      assemblies
      refdata

    main:
      makeChromosomeFastas(assemblies, refdata)

      // chrfasta=Channel.from("${makeChromosomeFastas.out.wd}/*.fasta")

      // println(chrfasta)

      // Channel.fromList(makeChromosomeFastas.out.chr_multifastas).view()

      chromosomeAlignment(makeChromosomeFastas.out.chr_multifastas)
}
