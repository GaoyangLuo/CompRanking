#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin

# get parent dir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# read YAML conda bin path
YAML_FILE="$PARENT_DIR/test_yaml.yaml"

# read YAML using python
CONDA_BIN_PATH=$(python3 -c "
import yaml
with open('$YAML_FILE', 'r') as f:
    data = yaml.safe_load(f)
print(data['CompRanking']['abs_path_to_conda_bin'])
")

# make sure get path
if [[ -z "$CONDA_BIN_PATH" ]]; then
    echo "Error: Failed to get Conda bin path from $YAML_FILE"
    exit 1
fi

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
source ${CONDA_BIN_PATH}/activate CompRanking_dvf_env
if [ -e checkdone/${PREFIX}.DVF.done ]; then
	echo "The first round phage prediction file existed..."
else
    #time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the first round phage prediction..."	
	#Running DVF 
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
    do
    python submodels/DeepVirFinder/dvf.py -i ${i} -o ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DVF -c ${THREADS} -l ${FITERLENGTH}
    done
	#finish Running DVF
	echo "[TIMESTAMP] $(date) Running the first round phage prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the first round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.DVF.done
fi
conda deactivate

# #run Seeker
# source ${CONDA_BIN_PATH}/activate CompRanking_seeker_env
# if [ -e ${PREFIX}.SEEKER.done ]; then
# 	echo "The second round phage prediction file existed..."
# else
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running the second round phage prediction..."	
# 	#Running Seeker
# 	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*.fa
#     do
#     predict-metagenome ${i}
#     done
# 	#finish Running Seeker
# 	echo "[TIMESTAMP] $(date) Running the second round phage prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running the second round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.SEEKER.done
# 	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
# fi


