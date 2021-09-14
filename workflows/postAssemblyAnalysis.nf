// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {makeChromosomeFastas} from '../modules/postAssemblyAnalysisModules.nf'

// define workflow
workflow postAssemblyAnalysis {

    take:
      assemblies
      refdata

    main:
      makeChromosomeFastas(assemblies, refdata)

}
