#!/usr/bin/env bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
CONDA_BIN_PATH=~/miniconda/bin
# FILTERLENTGH=0

while getopts "p:i:m:t:l:o" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		m) CONDA_BIN_PATH=${OPTARG};;
		t) THREADS=${OPTARG};; 
		l) FILTERLENTGH=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done

# echo ${OUTPUT_DIR} 
# echo ${INPUT_DIR}

#### STEP1 Ativate preprocessing env ####
source $CONDA_BIN_PATH/activate CompRanking_preprocessing_env
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_result
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/AMR/RGI
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/AMR/DeepARG
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/AMR/ARGranking
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DVF
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/DEF
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/plascad
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/Seeker
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/MGE/MobileOG
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/Virulence/VFDB
mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/Virulence/PATRIC

#mkdir -p ${INPUT_DIR}/preprocessing/5M-1K_contigs

#### Step2 Filtering contigs ####
#fiter 5M
if [ -e checkdone/${PREFIX}.5Mfilter.done ]; then
	echo "5M_filetered file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Filtering 5M contigs..."	
	#fiter 5M contigs
	for i in ${INPUT_DIR}/*fa
	do
	base=${i%%.f*}
	seqmagick convert --min-length ${FILTERLENTGH} ${i} ${base}_5M_contigs.fa
	done
	#finish 5M filtering
	echo "[TIMESTAMP] $(date) Filtering 5M contigs... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Filtering 5M contigs took $(($ENDTIME - $STARTTIME)) sec."
	mv ${INPUT_DIR}/*5M_contigs.fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs
    # cp ${INPUT_DIR}/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file
	
	touch checkdone/${PREFIX}.5Mfilter.done
fi

#filter 5M-1K
# if [ -e ${PREFIX}.5M-1Kfilter.done ]; then
# 	echo "5M-1K_filetered file existed..."
# else
# 	#time start
#     STARTTIME=$(date +%s)
# 	echo "[TIMESTAMP] $(date) Filtering 5M-1K contigs..."
# 	#filter 5M-1K congtigs
# 	for i in ${INPUT_DIR}/*5M*
# 	do
# 	base=${i%%5M_con*}
# 	seqmagick convert --max-length 1000 ${i} ${base}_5M-1K_contigs.fa
# 	done
# 	#finish 5M-1K filtering
# 	echo "[TIMESTAMP] $(date) Filtering 5M-1K contigs... Done"
# 	ENDTIME=$(date +%s)
# 	echo "[TIMER] Filtering 5M-1K contigs took $(($ENDTIME - $STARTTIME)) sec."
# 	touch ${PREFIX}.5M-1Kfilter.done
# fi
#mv filtered files
#mv ${INPUT_DIR}/*5M-1K_contigs.fa ${INPUT_DIR}/preprocessing/5M-1K_contigs

#### Step3 Run prodigal ####
if [ -e checkdone/${PREFIX}.ORF_prediction.done ]; then
	echo "faa file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Predicting ORFs with prodigal..."
	#predict ORFs with prodigal
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa
	do
	base=${i%%.f*}
	prodigal -i ${i} -o ${base}.gff -a ${base}.faa -d ${base}.fna -f gff -p meta -q
	# mkdir -p ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/split/${base} #ERR.contigs_5M_contigs
	done
	#finish ORFs prediction
	echo "[TIMESTAMP] $(date) Predicting ORFs with prodigal... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Predicting ORFs with prodigal took $(($ENDTIME - $STARTTIME)) sec."
	touch checkdone/${PREFIX}.ORF_prediction.done
fi
conda deactivate

#### Step4 Modify faa file ####
# cp ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*faa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file
# cp ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file


#### Step5 Building index ####
if [ -e checkdone/${PREFIX}.index_build.done ]; then
	echo "index file existed..."
else
	cp ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*faa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file
	cp ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/ori_file
	sed -i 's/[^>]*ID=//;s/;.*//;s/*//' ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*faa
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*gff
	do
	base=${i%%.gf*}
	sed -i '/^#/d' ${i}
	cut -f 1,9 ${i} |cut -d';' -f1| sed 's/ID=//' > ${base}.index
	done
	touch checkdone/${PREFIX}.index_build.done
fi