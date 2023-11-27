#!/bin/bash
#PBS -q workq
#PBS -N bwt
#PBS -l nodes=c001:ppn=32
module load miniconda3
workdir=/lomi_home/gaoyang/software/CompRanking/tmp_DSR/DSR/CompRanking_intermediate/preprocessing/5M_contigs/cov
workdir=`pwd`
cd $workdir

source activate  bowtie2

bowtie2-build ${1} ${1}_bwt

for i in *_1.f*
do
bowtie2 --very-sensitive --no-unal  -x ${1}_bwt -1 $i -2 ${i/_1/_2} -S ${i/_1.f*}.sam -p 32
done

conda deactivate

source activate samtools
for i in *sam
do
samtools view -bS $i > ${i/.sam/}.bam
samtools sort ${i/.sam/}.bam -o ${i/.sam/}.sorted.bam
done

conda deactivate

#对比对结果相似度进行筛选
source activate bamm
for i in *sorted.bam
do
bamm filter -b $i --percentage_id 0.99 --percentage_aln 0.9
done

conda deactivate

#统计每个样品的比对结果
source activate bbmap
for i in *filtered.bam
do
pileup.sh in=$i out=${i/.bam/}.cov rpkm=${i/.bam/}.rpkm overwrite=true
done

conda deactivate
