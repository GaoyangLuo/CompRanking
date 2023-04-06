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
#run DVF
source ${CONDA_BIN_PATH}/activate genomad
if [ -e checkdone/${PREFIX}.genomad.done ]; then
	echo "The MGE prediction file existed..."
else
    #time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the first round MGE prediction..."	
	#Running genomad
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
    do
    genomad end-to-end --cleanup --splits 8 ${i} ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/genomad_output /lomi_home/gaoyang/software/CompRanking/databases/genomad_db
    done
	#finish Running genomad
	echo "[TIMESTAMP] $(date) Running the first round phage prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the first round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.genomad.done
fi
conda deactivate

