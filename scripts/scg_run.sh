#!/bin/bash
set -e
set -m
#bash scripts/scg_run.sh -i xxx/CompRanking_intermediate/preprocessing/5M_contigs -p  -t 
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
source ${CONDA_BIN_PATH}/activate CompRanking_alignment_env


COV_DIR=${INPUT_DIR}/cov

#blastp scg
if [ -e checkdone/${PREFIX}.scg.done ]; then
	echo "scg identification file existed..."
else
	echo "Running scg identification..."
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running scg identification..."	
	#Running scg
	for i in ${INPUT_DIR}/*faa
	do
	base=${i/.faa/}
	#scg
    diamond blastp --query ${i} \
	       		   --db ${WORK_DIR}/submodels/scg_alignment/v1_scg/scg.seq.dmnd \
           		   --out ${base}_scg_Protein_dimond.txt \
                   --evalue 1e-10 \
                   --outfmt 6 \
                   --max-target-seqs 1 \
                   --threads ${THREADS}
	# diamond blastp --query ${i} --db ${WORK_DIR}/databases/SARG/SARG.db.fasta.dmnd --out ${base}_SARG_Protein_diamond.txt --evalue 1e-5 --outfmt 6 --threads ${THREADS} --max-target-seqs 1 
	done
	#finish Running VF prediction
	echo "[TIMESTAMP] $(date) Running scg identification... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running scg identification took $(($ENDTIME - $STARTTIME)) sec."
	# touch checkdone/${PREFIX}.ARG_ranking.done
	mv ${INPUT_DIR}/*_scg_Protein_dimond.txt ${COV_DIR}
fi
