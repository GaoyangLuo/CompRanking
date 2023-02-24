#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
WORK_DIR=$(dirname ${SCRIPT_PATH})
# echo ${SCRIPT_PATH}
# echo ${WORK_DIR}

while getopts "i:p:m:t:o" option; do
	case "${option}" in
		i) INPUT_DIR=${OPTARG};;
		p) PREFIX=${OPTARG};;
		m) CONDA_BIN_PATH=${OPTARG};;
		t) THREADS=${OPTARG};; 
        o) OUTPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done

#run plascad
source ${CONDA_BIN_PATH}/activate CompRanking_hmmer_env

if [ -e checkdone/${PREFIX}.PLASCAD.done ]; then
	echo "plascad file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running plascad prediction..."	
	#Running plascad
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*contigs.fa
    do
	echo ${i}
    python ${WORK_DIR}/submodels/plas_cad/plascad.py -i ${i}
	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*_id_* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad/
	rm -rf ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*out
    done
#     mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*_Conj_plasmids_id_out ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad/  
#     mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*_mob_unconj_plasmids_id_out ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad/
#     mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*_unmob_plasmids_id_out ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad/
    #   mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*_id_* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad/
	#finish Running plascad
	echo "[TIMESTAMP] $(date) Running plascad prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running plascad prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.PLASCAD.done
fi

rm -rf ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file/*out

