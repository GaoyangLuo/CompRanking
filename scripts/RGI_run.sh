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
RGI_DIR=$(dirname ${AMR_DIR_tmp})/AMR/RGI

#run RGI
source ${CONDA_BIN_PATH}/activate CompRanking_rgi_env
#input=/lomi_home/gaoyang/microplastic_test/metacompare_data/2_assembly/5M/faa
if [ -e checkdone/${PREFIX}.RGI.done ]; then
	echo "RGI file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running ARG prediction..."	
	#Running ARG prediction
	for i in ${INPUT_DIR}/*.faa
    do
    rgi main --input_sequence ${i} --output_file ${i%%.fa*}.RGI.out --local --clean --input_type protein --num_threads ${THREADS}
    done
	#finish Running ARG prediction
	echo "[TIMESTAMP] $(date) Running ARG prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running ARG prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.RGI.done
	mv ${INPUT_DIR}/*RGI.out* ${RGI_DIR}
fi


conda deactivate
