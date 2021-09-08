#!/bin/bash

container_list=("preprocessing" "parapipe")

for item in ${container_list[@]}; do
    sudo singularity build ${item}.sif singularity.${item}
done
