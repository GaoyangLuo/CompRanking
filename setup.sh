#!/bin/bash

#pathset
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})

#config conda evn
conda env create -f preprocessing.yaml
conda env create -f seeker.yaml
conda env create -f dvf.yaml
conda env create -f PlasFlow.yaml
conda env create -f rgi.yaml

#adjust seeker env
CONDA_BIN=`cat ${SCRIPT_PATH}/test_yaml.yaml |shyaml get-value CompRanking.abs_path_to_conda_bin`
CONDA_PATH=$(basename ${CONDA_BIN})
SEEKER_DIR=${CONDA_PATH}/envs/seeker/lib/python3.7/site-packages/seeker
rm ${SEEKER_DIR}/command_line.py && cp scripts/command_line.py ${SEEKER_DIR}

#check 


