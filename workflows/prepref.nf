// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {getRefData} from '../modules/preprefModules.nf'
include {indexRefData} from '../modules/preprefModules.nf'

// define workflow
workflow prepRef {

    take:
      ref_id

    main:
      getRefData(ref_id)

      indexRefData(getRefData.out.refdata, ref_id)

    emit:
      refdata = getRefData.out.refdata
      ref_bt2index = indexRefData.out.bt2_index
      ref_fai = indexRefData.out.reffaidx
}
