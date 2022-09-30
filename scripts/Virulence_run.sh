#!/bin/bash
set -e
set -m

# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin
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
VIR_DIR_tmp=$(dirname ${INPUT_DIR})
VIR_DIR=$(dirname ${VIR_DIR_tmp})/Virulence
echo ${CONDA_BIN_PATH}
echo ${INPUT_DIR}
###############################Virulence Factore Prediction######################
#run HMMersearch
#activate environment
source ${CONDA_BIN_PATH}/activate hmmer
#hmmer
if [ -e ${PREFIX}.hmm.done ]; then
	echo "hmm predition file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running hmm prediction..."	
	#Running Virulence Factor prediction
	for i in ${INPUT_DIR}/*faa
	do
	base=${i%%.f*}
	hmmsearch --cpu ${THREADS} --noali --notextw --tblout ${base}.hmmscan databases/virulence/Virulence_factor.hmm ${i} 
	done
	#finish Running VF prediction
	echo "[TIMESTAMP] $(date) Running hmm prediction... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running hmm prediction took $(($ENDTIME - $STARTTIME)) sec."
	touch ${PREFIX}.hmm.done
fi 

# mv ${INPUT_DIR}/*hmmscan ${VIR_DIR}

#VFDB
# source ${CONDA_BIN_PATH}/activate argranker

# if [ -e ${PREFIX}.VFDB.done ]; then
# 	echo "Virulence predition file existed..."
# else
# 	echo "Running VFDB prediction..."
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running VFDB prediction..."	
# 	#Running Virulence Factor prediction
# 	for i in ${INPUT_DIR}/*faa
# 	do
# 	base=${i%%.f*}
# 	diamond blastp --query ${i} --db ${WORK_DIR}/databases/VFDB/VFDB_setB_pro.fas.dmnd --out ${base}_VFDB.out --evalue 1e-10 --outfmt 6 --threads ${THREADS} 
# 	done
# 	#finish Running VF prediction
# 	echo "[TIMESTAMP] $(date) Running VFDB prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running VFDB prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.VFDB.done
# fi
# conda deactivate
# mv ${INPUT_DIR}/*out ${VIR_DIR}

#PATRIC
# source ${CONDA_BIN_PATH}/activate argranker

# if [ -e ${PREFIX}.PATRIC.done ]; then
# 	echo "Virulence predition file existed..."
# else
# 	echo "Running PATRIC prediction..."
# 	#time start
# 	STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Running PATRIC prediction..."	
# 	#Running Virulence Factor prediction
# 	for i in ${INPUT_DIR}/*fa
# 	do
# 	base=${i%%.f*}
# 	blastn -query ${i} -db ${WORK_DIR}/databases/PATRIC/PATRIC -out ${base}_PATRIC.out -evalue 1e-10 -outfmt 6 -num_threads ${THREADS}
# 	done
# 	#finish Running VF prediction
# 	echo "[TIMESTAMP] $(date) Running PATRIC prediction... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Running PATRIC prediction took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.PATRIC.done
# fi
# mv ${INPUT_DIR}/*PATRIC.out ${VIR_DIR}
# hmmsearch --cpu 16 --noali --notextw --tblout ERR1191817_5M.hmmscan /lomi_home/gaoyang/software/CompRanking/databases/virulence/Virulence_factor.hmm /lomi_home/gaoyang/software/CompRanking/test/CompRanking/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.faa
# hmmsearch --cpu 16 --noali --notextw --tblout ERR1191817_5M.hmmscan /lomi_home/gaoyang/software/CompRanking/databases/virulence/Virulence_factor.hmm /lomi_home/gaoyang/software/CompRanking/test/CompRanking/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.faa
# diamond blastp --query /lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817.contigs_5M_contigs.faa --db /lomi_home/gaoyang/software/CompRanking/databases/VFDB/VFDB_setA_pro.fas.dmnd --out /lomi_home/gaoyang/software/CompRanking/test/CompRanking/CompRanking_intermediate/preprocessing/5M_contigs/ERR1191817_VFDB.txt --evalue 1e-5 --outfmt 6 --max-target-seqs 1 --threads 20