#!/bin/bash
#PBS -q workq
#PBS -N bwt
#PBS -l nodes=c001:ppn=32
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
    if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa`" = "" ]; then 
        echo "No fasta files"
        ln -s ${INPUT_DIR}/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    else 
        echo "clean old fasta record..."
        rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
    fi

    if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq`" = "" ]; then 
        echo "No fastq files"
        ln -s ${INPUT_DIR}/*fq ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    else 
        echo "clean old fastq record..."
        rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq
    fi

    if [ "`ls -A ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna`" = "" ]; then 
        echo "No fna files"
        ln -s ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fna ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    else 
        echo "clean old fna record..."
        rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna
    fi
    

    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fq
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
    # rm ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fna
    
    #ln fastq fasta
    # ln -s ${INPUT_DIR}/*fq ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # ln -s ${INPUT_DIR}/*fa ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov
    # ln -s ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/*fna ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov

    #write a for circle to run all the files in the cov dir
	for i in ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov/*fa
    do
    base=${i/.fa/}

    #run bowtie2
    bowtie2-build ${base}_5M_contigs.fna ${base}_5M_contigs.fna_bwt
    bowtie2 --very-sensitive \
            --no-unal \
            -x ${base}_5M_contigs.fna_bwt \
            -1 ${base}_1.fq -2 ${base}_2.fq \
            -S ${base}.sam \
            -p ${THREADS}

    #run samtools
    samtools view -bS ${base}.sam > ${base}.bam
    samtools sort ${base}.bam -o ${base}.sorted.bam
    
    #run bamm
    echo "[TIMESTAMP] $(date) Running bamm..."
    bamm filter -b ${base}.sorted.bam -o ${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs/cov --percentage_id 0.99 --percentage_aln 0.9

    #run bbmap
    echo "[TIMESTAMP] $(date) Running bbmap..."
    pileup.sh in=${base}.sorted_filtered.bam out=${base}_5M_contigs_gene.cov rpkm=${base}_5M_contigs_gene.rpkm overwrite=true

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
# source ~/miniconda/bin/activate CompRanking_abundance_env
# workdir=/lomi_home/gaoyang/software/CompRanking/tmp_NPSW/NPSW/CompRanking_intermediate/preprocessing/5M_contigs/cov
# workdir=`pwd`
# cd $workdir

# bowtie2-build bowtie2-build ${1} ${1}_bwt

# for i in *_1.f*
# do
# bowtie2 --very-sensitive --no-unal  -x ${1}_bwt -1 $i -2 ${i/_1/_2} -S ${i/_1.f*}.sam -p 32
# done

# for i in *sam
# do
# samtools view -bS $i > ${i/.sam/}.bam
# samtools sort ${i/.sam/}.bam -o ${i/.sam/}.sorted.bam
# done

# for i in *sorted.bam
# do
# bamm filter -b $i --percentage_id 0.99 --percentage_aln 0.9
# done

# for i in *filtered.bam
# do
# pileup.sh in=$i out=${i/.bam/}.cov rpkm=${i/.bam/}.rpkm overwrite=true
# done



#################################################
# module load miniconda3
# workdir=/lomi_home/gaoyang/software/CompRanking/tmp_NPSW/NPSW/CompRanking_intermediate/preprocessing/5M_contigs/cov
# workdir=`pwd`
# cd $workdir

# source activate bowtie2

# bowtie2-build ${1} ${1}_bwt

# for i in *_1.f*
# do
# bowtie2 --very-sensitive --no-unal  -x ${1}_bwt -1 $i -2 ${i/_1/_2} -S ${i/_1.f*}.sam -p 32
# done

# conda deactivate

# source activate samtools
# for i in *sam
# do
# samtools view -bS $i > ${i/.sam/}.bam
# samtools sort ${i/.sam/}.bam -o ${i/.sam/}.sorted.bam
# done

# conda deactivate

# #对比对结果相似度进行筛选
# source activate bamm
# for i in *sorted.bam
# do
# bamm filter -b $i --percentage_id 0.99 --percentage_aln 0.9
# done

# conda deactivate

# #统计每个样品的比对结果
# source activate bbmap
# for i in *filtered.bam
# do
# pileup.sh in=$i out=${i/.bam/}.cov rpkm=${i/.bam/}.rpkm overwrite=true
# done

# conda deactivate
