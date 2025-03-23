#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
CONDA_BIN_PATH=~/miniconda/bin
# FITERLENGTH=500

while getopts "p:i:m:t:l:o" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		m) CONDA_BIN_PATH=${OPTARG};;
		t) THREADS=${OPTARG};; 
		l) FITERLENGTH=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done



#run DEF
source ${CONDA_BIN_PATH}/activate def
if [ -e checkdone/${PREFIX}.DEF.done ]; then
	echo "The DEF file existed..."
else
    #time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the DEF prediction..."	
	#Running DVF 
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
    do
    python submodels/DeepMicrobeFinder/predict.py -i ${i} -e one-hot -d submodels/DeepMicrobeFinder/models/one-hot-models/ -m hybrid -o ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DEF 
    done
	#finish Running DVF
	echo "[TIMESTAMP] $(date) Running the DEF prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the DEF prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.DEF.done
fi
conda deactivate



