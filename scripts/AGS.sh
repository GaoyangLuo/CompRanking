#!/bin/bash
#PBS -N MicrobeCensus
#PBS -l nodes=c005:ppn=32


threads=32
#export PATH=$PATH:~/software/MicrobeCensus/scripts
source /lomi_home/gaoyang/miniconda/bin/activate base
#export PYTHONPATH=$PYTHONPATH:~/software/MicrobeCensus

workdir=/lomi_home/gaoyang/software/CompRanking/tmp_NPSW

cd $workdir
for i in *_1.fq
do
run_microbe_census.py -t $threads ${i},${i/_1.fq/_2.fq} ${i/_1.fq}.contigs.AGS.txt
done

conda deactivate