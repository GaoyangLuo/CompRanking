#!/bin/bash

# KKDB=/lomi_home/gaoyang/db/kraken2/202203
# workdir=/lomi_home/gaoyang/microplastic_test/dashariver/metagenome_illunina
# cd $workdir

# threads=32

# #create_resipotory_to_store_rawdata
# mkdir arg_run_flexbar
# input=$workdir/arg_run_flexbar
# #linkdata
# ln -s /lomi_home/gaoyang/microplastic_test/dashariver/metagenome_illunina/work/clean_reads/*fq ${input}
# #activate argranker
# source /lomi_home/gaoyang/miniconda/bin/activate argranker
# #follow_steps1
# arg_ranker -i ${input} -t 32 -kkdb /lomi_home/gaoyang/db/kraken2/202203
# #follow_steps2
# sh $workdir/arg_ranking/script_output//arg_ranker.sh
set -e
set -m

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
source ${CONDA_BIN_PATH}/activate CompRanking_alignment_env

# if [ -e ${PREFIX}.ARG_ranking.done ]; then
# 	echo "ARG_ranking predition file existed..."
# else
# 	echo "Running ARG_ranking prediction..."
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running ARG_ranking prediction..."	
# 	#Running Virulence Factor prediction
# 	for i in ${INPUT_DIR}/*fa
# 	do
# 	base=${i%%.f*}
# 	diamond blastx --query ${i} --db ${WORK_DIR}/databases/SARG/SARG.db.fasta.dmnd --out ${base}_SARG_diamond.txt --evalue 1e-5 --outfmt 6 --threads ${THREADS} --max-target-seqs 1
# 	done
# 	#finish Running VF prediction
# 	echo "[TIMESTAMP] $(date) Running ARG_ranking prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running ARG_ranking prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.ARG_ranking.done
# fi

AMR_DIR_tmp=$(dirname ${INPUT_DIR})
SARG_DIR=$(dirname ${AMR_DIR_tmp})/AMR/ARGranking
MobileOG_DIR=$(dirname ${AMR_DIR_tmp})/MGE/MobileOG

#blastp SARG
if [ -e checkdone/${PREFIX}.ARG_ranking.done ]; then
	echo "ARG_ranking predition file existed..."
else
	echo "Running ARG_ranking prediction..."
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running ARG_ranking prediction..."	
	#Running Virulence Factor prediction
	for i in ${INPUT_DIR}/*faa
	do
	base=${i%%.f*}
	diamond blastp --query ${i} --db ${WORK_DIR}/databases/SARG/SARG.db.fasta.dmnd --out ${base}_SARG_Protein_diamond.txt --evalue 1e-5 --outfmt 6 --threads ${THREADS} --max-target-seqs 1 
	done
	#finish Running VF prediction
	echo "[TIMESTAMP] $(date) Running ARG_ranking prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running ARG_ranking prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.ARG_ranking.done
	mv ${INPUT_DIR}/*_SARG_Protein_diamond.txt ${SARG_DIR}
fi


#blasp 