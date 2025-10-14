#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin
# get the parent directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# read  conda bin path from test_yaml.yaml
YAML_FILE="$PARENT_DIR/test_yaml.yaml"

#  use Python read YAML
CONDA_BIN_PATH=$(python3 -c "
import yaml
with open('$YAML_FILE', 'r') as f:
    data = yaml.safe_load(f)
print(data['CompRanking']['abs_path_to_conda_bin'])
")

# Check if CONDA_BIN_PATH is empty
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



#run DEF
source ${CONDA_BIN_PATH}/activate CompRanking_def_env
if [ -e checkdone/${PREFIX}.DEF.done ]; then
	echo "The DEF file existed..."
else
    #time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the DEF prediction..."	
	#Running DVF 
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
    do
    DeepMicroClass predict -i ${i} -e onehot -md hybrid --cpu_thread $THREADS -o ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DEF -d cpu 
    done
	#finish Running DVF
	echo "[TIMESTAMP] $(date) Running the DEF prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the DEF prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.DEF.done
fi
conda deactivate



