#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
CONDA_BIN_PATH=~/miniconda/bin

while getopts "p:i:m:t:o" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		m) CONDA_BIN_PATH=${OPTARG};;
		t) THREADS=${OPTARG};; 
        o) OUTPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done


#run Seeker
source ${CONDA_BIN_PATH}/activate CompRanking_seeker_env
if [ -e checkdone/${PREFIX}.SEEKER.done ]; then
	echo "The second round phage prediction file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the second round phage prediction..."	
	#Running Seeker
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*.fa
    do
    predict-metagenome ${i}
    done
	#finish Running Seeker
	echo "[TIMESTAMP] $(date) Running the second round phage prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the second round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.SEEKER.done
	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
fi