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

#run plasflow
# source ${CONDA_BIN_PATH}/activate plasflow

# if [ -e ${PREFIX}.PLASFLOW.done ]; then
# 	echo "plasflow file existed..."
# else
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running plasmid prediction..."	
# 	#Running plasflow
# 	for i in ${INPUT_DIR}/${PREFIX}/preprocessing/5M_contigs/*fa
#     do
#     base=${i%%.f*}
#     PlasFlow.py --input ${i} --output ${base}_plasflow_predictions.tsv --threshold 0.7
#     done
# 	#finish Running plasflow
# 	echo "[TIMESTAMP] $(date) Running plasmid prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running plasmid prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.PLASFLOW.done
# 	mv ${INPUT_DIR}/${PREFIX}/preprocessing/5M_contigs/*plasflow* ${INPUT_DIR}/${PREFIX}/MGE/Plasflow
# fi
# conda deactivate


#run DVF
source ${CONDA_BIN_PATH}/activate dvf
if [ -e ${PREFIX}.DVF.done ]; then
	echo "The first round phage prediction file existed..."
else
    #time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the first round phage prediction..."	
	#Running DVF 
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
    do
    python submodels/DeepVirFinder/dvf.py -i ${i} -o ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DVF -c 16 -l 500
    done
	#finish Running DVF
	echo "[TIMESTAMP] $(date) Running the first round phage prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the first round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch ${PREFIX}.DVF.done
fi
conda deactivate

#run Seeker
source ${CONDA_BIN_PATH}/activate seeker
if [ -e ${PREFIX}.SEEKER.done ]; then
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
	touch ${PREFIX}.SEEKER.done
	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
fi


