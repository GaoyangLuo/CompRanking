#!/bin/bash

#config conda evn
conda env create -f preprocessing.yaml
conda env create -f seeker.yaml
conda env create -f dvf.yaml
conda env create -f PlasFlow.yaml
conda env create -f rgi.yaml

#adjust seeker env
CONDA_BIN=/lomi_home/gaoyang/miniconda/bin
CONDA_PATH=$(basename ${CONDA_BIN})
SEEKER_DIR=${CONDA_PATH}/CompRanking_seeker_env/lib/python3.7/site-packages
rm ${SEEKER_DIR}/seeker/command_line.py && cp scripts/command_line.py ${SEEKER_DIR}/seeker/

#check 


