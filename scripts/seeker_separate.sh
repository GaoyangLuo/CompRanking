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


# Step1 创建有名管道
[ -e ./fd1 ] || mkfifo ./fd1

# 创建文件描述符，以可读（<）可写（>）的方式关联管道文件，这时候文件描述符3就有了有名管道文件的所有特性
exec 3<> ./fd1   

# 关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
rm -rf ./fd1                    

# Step2 创建令牌 
for i in `seq 1 12`;
do
    # echo 每次输出一个换行符,也就是一个令牌
    echo >&3                   
done


#run Seeker
source ${CONDA_BIN_PATH}/activate CompRanking_seeker_env
if [ -e checkdone/${PREFIX}.SEEKER.done ]; then
	echo "The second round phage prediction file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running the second round phage prediction..."	
	#Running Seeker
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*.fa
	do
	
	read -u3                    # read 命令每次读取一行，也就是拿到一个令牌 
	{
		sleep 5
		echo "${i} 正在执行phage预测"
		predict-metagenome ${i}
		echo >&3                # 执行完一条命令会将令牌放回管道
		
	}&
	done
	wait

	exec 3<&-                       # 关闭文件描述符的读
	exec 3>&-                       # 关闭文件描述符的写

	#finish Running Seeker
	echo "[TIMESTAMP] $(date) Running the second round phage prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running the second round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.SEEKER.done
	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
fi




# 	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*.fa
#     do
#     predict-metagenome ${i}
#     done
# 	#finish Running Seeker
# 	echo "[TIMESTAMP] $(date) Running the second round phage prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running the second round phage prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch checkdone/${PREFIX}.SEEKER.done
# 	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
# fi





#run Seeker
source ${CONDA_BIN_PATH}/activate CompRanking_seeker_env
if [ -e checkdone/${PREFIX}.SEEKER.done ]; then
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
	touch checkdone/${PREFIX}.SEEKER.done
	mv ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/seeker* ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
fi