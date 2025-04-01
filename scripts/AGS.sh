#!/bin/bash
#SBATCH -p jufeng   
#SBATCH -q huge                       
#SBATCH -J microCensus              
#SBATCH -c 16


threads=16
#export PATH=$PATH:~/software/MicrobeCensus/scripts

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

# source $CONDA_BIN_PATH/activate CompRanking
#export PYTHONPATH=$PYTHONPATH:~/software/MicrobeCensus


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


source $CONDA_BIN_PATH/activate microbecensus_py2

workdir=$INPUT_DIR

cd $workdir
for i in *_1.fq
do
run_microbe_census.py -t $threads ${i},${i/_1.fq/_2.fq} ${i/_clean_1.fq}.contigs.AGS.txt
done

mv $workdir/*AGS* $workdir/$PREFIX/CompRanking_intermediate/preprocessing/5M_contigs/AGS

conda deactivate