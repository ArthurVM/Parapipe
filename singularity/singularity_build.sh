#!/bin/bash

container_list=("preprocessing" "parapipe" "moi")

for item in ${container_list[@]}; do
    sudo singularity build ${item}.sif singularity.${item}
done
