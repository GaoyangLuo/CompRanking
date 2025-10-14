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
source ${CONDA_BIN_PATH}/activate CompRanking_deeparg_env
#input=/lomi_home/gaoyang/microplastic_test/metacompare_data/2_assembly/5M/faa
if [ -e checkdone/${PREFIX}.DeepARG.done ]; then
	echo "DeepARG file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running DeepARG prediction..."	
	#Running ARG prediction
	for i in ${INPUT_DIR}/*.faa
    do
    deeparg predict --model LS --type prot --input ${i} --out ${i%%.fa*}_DeepARG.out -d databases/deepargdata1.0.2 --arg-alignment-overlap 0.7
    done
	#finish Running ARG prediction
	echo "[TIMESTAMP] $(date) Running DeepARG prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running ARG prediction took $(($ENDTIME - $STARTTIME)) sec."
	mv ${INPUT_DIR}/*DeepARG.out* ${DeepARG_DIR}
	touch checkdone/${PREFIX}.DeepARG.done
fi



conda deactivate





