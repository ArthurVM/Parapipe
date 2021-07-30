// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {DownloadRefData} from '../modules/preprefModules.nf'
include {IndexRefData} from '../modules/preprefModules.nf'

// define workflow
workflow prepRef {

    take:
      genome_id

    main:
    DownloadRefData(genome_id)
    IndexRefData(DownloadRefData.out.refdata, genome_id)

    emit:
    ref_bt2index = IndexRefData.out.indexBT2
    ref_fai = IndexRefData.out.reffaidx
    ref_dict = IndexRefData.out.refdict
}
