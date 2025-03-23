#!/bin/bash
set -e
set -m

#######This scripts is for ARG abundance benchmark################
#######by using DeepARG to predict ARGs from ORFs#################
#######prediced from the dna sequences files using prodigal#######

# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
WORK_DIR=$(dirname ${SCRIPT_PATH})

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
mkdir -p $(dirname ${AMR_DIR_tmp})/AMR/DeepARG/benchmark
DeepARG_DIR=$(dirname ${AMR_DIR_tmp})/AMR/DeepARG/benchmark

#run DeepARG
source ${CONDA_BIN_PATH}/activate CompRanking_deeparg_env
#input=preprocessing/5M_contigs/fna
rm checkdone/${PREFIX}.DeepARG.done
if [ -e checkdone/${PREFIX}.DeepARG.done ]; then
	echo "DeepARG file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running DeepARG prediction..."	
	#Running ARG prediction
	for i in ${INPUT_DIR}/*.fna
    do
    deeparg predict --model LS --type nucl --input ${i} --out ${input}/${i%%.fa*}_DeepARG.out -d databases/deepargdata1.0.2 --arg-alignment-overlap 0.7
    done
	#finish Running ARG prediction
	echo "[TIMESTAMP] $(date) Running DeepARG prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running ARG prediction took $(($ENDTIME - $STARTTIME)) sec."
	mv ${INPUT_DIR}/*DeepARG.out* ${DeepARG_DIR}
	touch checkdone/${PREFIX}.DeepARG.done
fi



conda deactivate





