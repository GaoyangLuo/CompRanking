#!/bin/bash

#pathset
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})

#config conda evn
conda env create -f envs/alignment.yaml
conda env create -f envs/deeparg.yaml
conda env create -f envs/preprocessing.yaml
conda env create -f envs/seeker.yaml
conda env create -f envs/dvf.yaml
conda env create -f envs/rgi.yaml
conda env create -f envs/hmmer.yaml
conda env create -f envs/def.yaml

#adjust seeker env
CONDA_BIN=`cat ${SCRIPT_PATH}/test_yaml.yaml |shyaml get-value CompRanking.abs_path_to_conda_bin`
CONDA_PATH=$(dirname ${CONDA_BIN})
SEEKER_DIR=${CONDA_PATH}/envs/CompRanking_seeker_env/lib/python3.7/site-packages/seeker
rm ${SEEKER_DIR}/command_line.py && cp ${SCRIPT_PATH}/scripts/command_line.py ${SEEKER_DIR}

#make checkdone
mkdir checkdone



