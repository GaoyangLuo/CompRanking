#!/bin/bash
#PBS -q workq
#PBS -N bwt
#PBS -l nodes=c004:ppn=22
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
CONDA_BIN_PATH=~/miniconda/bin
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
WORK_DIR=$(dirname ${SCRIPT_PATH})

while getopts "p:i:m:t:o:r" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		i) INPUT_DIR=${OPTARG};;
		m) CONDA_BIN_PATH=${OPTARG};;
		t) THREADS=${OPTARG};; 
        o) OUTPUT_DIR=${OPTARG};;
		*) exit 1;;
	esac
done

#run cov
source ${CONDA_BIN_PATH}/activate CompRanking_abundance_env
#input=/lomi_home/gaoyang/microplastic_test/metacompare_data/2_assembly/5M/faa
if [ -e checkdone/${PREFIX}.cov.done ]; then
	echo "cov file existed..."
else
	#time start
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running cov..."	
	#Running cov prediction
    # if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa`" = "" ]; then 
    #     echo "No fasta files"
    #     ln -s ${INPUT_DIR}/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # else 
    #     echo "clean old fasta record..."
    #     rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
    # fi

    # if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq`" = "" ]; then 
    #     echo "No fastq files"
    #     ln -s ${INPUT_DIR}/*fq ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # else 
    #     echo "clean old fastq record..."
    #     rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq
    #     ln -s ${INPUT_DIR}/*fq ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # fi

    # if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna`" = "" ]; then 
    #     echo "No fna files"
    #     ln -s ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fna ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # else 
    #     echo "clean old fna record..."
    #     rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna
    #     ln -s ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fna ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # fi
    

    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna
    
    #ln fastq fasta
    # ln -s ${INPUT_DIR}/*fq ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # ln -s ${INPUT_DIR}/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # ln -s ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fna ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov

    #write a for circle to run all the files in the cov dir
    # bowtie2-build ${INPUT_DIR}/crassphage.fasta ${INPUT_DIR}/crassphage.fasta_bwt

	for i in ${INPUT_DIR}/*fa
    do
    base=${i/.fa/}

    #run bowtie2
    bowtie2 --very-sensitive \
            --no-unal \
            -x ${INPUT_DIR}/crassphage.fasta_bwt \
            -1 ${base}_1.fq -2 ${base}_2.fq \
            -S ${base}_crassphage.sam \
            -p ${THREADS}

    #run samtools
    samtools view -bS ${base}_crassphage.sam > ${base}_crassphage.bam
    samtools sort ${base}_crassphage.bam -o ${base}_crassphage.sorted.bam
    
    #run bamm
    echo "[TIMESTAMP] $(date) Running bamm..."
    bamm filter -b ${base}_crassphage.sorted.bam -o ${INPUT_DIR} --percentage_id 0.95 --percentage_aln 0.95

    #run bbmap
    echo "[TIMESTAMP] $(date) Running bbmap..."
    pileup.sh in=${base}_crassphage.sorted_filtered.bam out=${base}_5M_contigs_crassphage.cov rpkm=${base}_5M_contigs_crassphage.rpkm overwrite=true

    done
	#finish Running cov prediction
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
	echo "[TIMESTAMP] $(date) Running cov... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running cov took $(($ENDTIME - $STARTTIME)) sec."
	# mv ${INPUT_DIR}/*DeepARG.out* ${DeepARG_DIR}
	# touch checkdone/${PREFIX}.cov.done
fi

##################################
