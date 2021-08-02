// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {downloadRefData} from '../modules/preprefModules.nf'
include {indexRefData} from '../modules/preprefModules.nf'

// define workflow
workflow prepRef {

    take:
      genome_id

    main:
    downloadRefData(genome_id)
    indexRefData(downloadRefData.out.refdata, genome_id)

    emit:
    ref_bt2index = indexRefData.out.bt2_index
    ref_fai = indexRefData.out.reffaidx
    ref_dict = indexRefData.out.refdict
}
