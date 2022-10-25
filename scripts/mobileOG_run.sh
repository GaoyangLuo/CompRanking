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
CONDA_BIN_PATH=~/miniconda/bin
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
source ${CONDA_BIN_PATH}/activate argranker

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

#blastp mobileOG
if [ -e ${PREFIX}.MobileOG.done ]; then
	echo "MobileOG predition file existed..."
else
	echo "Running MobileOG prediction..."
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running MobileOG prediction..."	
	#Running Virulence Factor prediction
	for i in ${INPUT_DIR}/*faa
	do
	base=${i%%.f*}
	diamond blastp --query ${i} --db ${WORK_DIR}/databases/MobileOG-db/mobileOG-db_aa.dmnd --out ${base}_mobileOG_diamond.txt --evalue 1e-10 --outfmt 6 --threads ${THREADS} --max-target-seqs 1
	done
	#finish Running VF prediction
	echo "[TIMESTAMP] $(date) Running MobileOG prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running MobileOG prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch ${PREFIX}.MobileOG.done
	mv ${INPUT_DIR}/*mobileOG_diamond.txt ${MobileOG_DIR}
fi

