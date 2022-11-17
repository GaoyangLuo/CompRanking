#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
CONDA_BIN_PATH=~/miniconda/bin
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
WORK_DIR=$(dirname ${SCRIPT_PATH})

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

#######################################Preset####################################
AMR_DIR_tmp=$(dirname ${INPUT_DIR})
DeepARG_DIR=$(dirname ${AMR_DIR_tmp})/AMR/DeepARG

#run DeepARG
source ${CONDA_BIN_PATH}/activate deeparg1.0.2
#input=/lomi_home/gaoyang/microplastic_test/metacompare_data/2_assembly/5M/faa
if [ -e ${PREFIX}.DeepARG.done ]; then
	echo "DeepARG file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running DeepARG prediction..."	
	#Running ARG prediction
	for i in ${INPUT_DIR}/*.faa
    do
    deeparg predict --model LS --type prot --input ${i} --out ${input}/${i%%.fa*}_DeepARG.out -d databases/deepargdata1.0.2 --arg-alignment-overlap 0.7
    done
	#finish Running ARG prediction
	echo "[TIMESTAMP] $(date) Running DeepARG prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running ARG prediction took $(($ENDTIME - $STARTTIME)) sec."
	mv ${INPUT_DIR}/*DeepARG.out* ${DeepARG_DIR}
	touch ${PREFIX}.DeepARG.done
fi



conda deactivate





