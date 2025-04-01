#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin
# 获取当前脚本所在目录的上一级目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# 读取 YAML 文件中的 conda bin 路径
YAML_FILE="$PARENT_DIR/test_yaml.yaml"

# 使用 Python 解析 YAML
CONDA_BIN_PATH=$(python3 -c "
import yaml
with open('$YAML_FILE', 'r') as f:
    data = yaml.safe_load(f)
print(data['CompRanking']['abs_path_to_conda_bin'])
")

# 确保获取到路径
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



