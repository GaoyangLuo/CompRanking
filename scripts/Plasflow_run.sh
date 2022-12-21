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
# source ${CONDA_BIN_PATH}/activate CompRanking_plasflow_env #plasflow

# if [ -e ${PREFIX}.PLASFLOW.done ]; then
# 	echo "plasflow file existed..."
# else
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running plasmid prediction..."	
# 	#Running plasflow
# 	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
#     do
#     base=${i%%.f*}
#     PlasFlow.py --input ${i} --output ${base}_plasflow_predictions.tsv --threshold 0.7
#     done
# 	#finish Running plasflow
# 	echo "[TIMESTAMP] $(date) Running plasmid prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running plasmid prediction took $(($ENDTIME - $STARTTIME)) sec."
	
# 	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*plasflow* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow
# 	touch ${PREFIX}.PLASFLOW.done
# fi
# conda deactivate

#run plasflow2
source ${CONDA_BIN_PATH}/activate CompRanking_plasflow_env #plasflow

if [ -e checkdone/${PREFIX}.PLASFLOW.done ]; then
	echo "plasflow file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running plasmid prediction..."	
	#Running plasflow
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/split/*contigs
	do
	base_name=$(basename $i)
	echo $i
	echo $base_name
		for j in ${i}/*fasta
    	do
		base=${j%%.fas*}
    	PlasFlow.py --input ${j} --output ${base}_plasflow_predictions.tsv --threshold 0.7
		mv ${base}_plasflow_predictions.tsv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow
		rm ${i}/*plasflow*fasta
    	done
	cat ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/*_plasflow_predictions.tsv > ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/tmp.tsv
	cut -f 3,6 ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/tmp.tsv > ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/${base_name}_plasflow_predictions_final.tsv
	rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/*_plasflow_predictions.tsv
	rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow/tmp.tsv
	done
	#finish Running plasflow
	echo "[TIMESTAMP] $(date) Running plasmid prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running plasmid prediction took $(($ENDTIME - $STARTTIME)) sec."
	
	# mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*plasflow* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Plasflow
	touch checkdone/${PREFIX}.PLASFLOW.done
fi
conda deactivate



# #run DVF
# source ~/miniconda/bin/activate dvf
# if [ -e ${PREFIX}.DVF.done ]; then
# 	echo "The first round phage prediction file existed..."
# else
#     #time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running the first round phage prediction..."	
# 	#Running DVF 
# 	for i in ${INPUT_DIR}/${PREFIX}/preprocessing/5M_contigs/*.fa
#     do
#     python submodels/DeepVirFinder/dvf.py -i ${i} -o ${INPUT_DIR}/CompRanking/MGE/DVF -c 16 -l 500
#     done
# 	#finish Running DVF
# 	echo "[TIMESTAMP] $(date) Running the first round phage prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running the first round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.DVF.done
# fi
# conda deactivate

# #run Seeker
# source ${CONDA_BIN_PATH}/activate seeker
# if [ -e ${PREFIX}.SEEKER.done ]; then
# 	echo "The second round phage prediction file existed..."
# else
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running the second round phage prediction..."	
# 	#Running Seeker
# 	for i in ${INPUT_DIR}/${PREFIX}/preprocessing/5M_contigs/*.fa
#     do
#     predict-metagenome ${i}
#     done
# 	#finish Running Seeker
# 	echo "[TIMESTAMP] $(date) Running the second round phage prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running the second round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.SEEKER.done
# fi
# mv ${INPUT_DIR}/${PREFIX}/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/MGE/Seeker

