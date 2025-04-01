#!/bin/bash
#SBATCH -p jufeng   
#SBATCH -q huge                       
#SBATCH -J prodigal           
#SBATCH -c 64

# source ~/miniconda3/bin/activate vs2
# mamba init
# mamba activate vs2
source ~/miniconda3/bin/activate CompRanking_preprocessing_env

WORKDIR=/storage/jufengLab/luogaoyang/software/CompRanking-main/tmp_soil_distinct/distinct/CompRanking_intermediate/preprocessing/5M_contigs/fna2faa_

cd $WORKDIR

# Step1 创建有名管道
[ -e ./fd1 ] || mkfifo ./fd1

# 创建文件描述符，以可读（<）可写（>）的方式关联管道文件，这时候文件描述符3就有了有名管道文件的所有特性 
exec 3<> ./fd1   

# 关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
rm -rf ./fd1                    

# Step2 创建令牌 
for i in `seq 1 10`;
do
    # echo 每次输出一个换行符,也就是一个令牌
    echo >&3                   
done

#run multi-threads hmm
STARTTIME=$(date +%s)
echo "[TIMESTAMP] $(date) Running the gene second round phage prediction..."	
#Running Seeker
for i in *fna
do
filename=${i%%.fna*}

read -u3                    # read 命令每次读取一行，也就是拿到一个令牌 
{
  sleep 3
  echo "${i} 正在执行prodigal预测"

  prodigal -i ${i} -a ${filename}.fna2faa.faa  -p meta -q

  echo >&3                # 执行完一条命令会将令牌放回管道
  
}&
done
wait

exec 3<&-                       # 关闭文件描述符的读
exec 3>&-                       # 关闭文件描述符的写

#finish Running Seeker
echo "[TIMESTAMP] $(date) progdigal运行... Done"
ENDTIME=$(date +%s)
echo "[TIMER] Running the progdigal took $(($ENDTIME - $STARTTIME)) sec."


# for i in /storage/jufengLab/luogaoyang/metagenome_project/DSR/2_assembly/*fa
# do
# base=${i%%.con*}
# prodigal -i ${i} -d ${base}.fna -p meta
# done

# prodigal -i ${i} -o ${filebase}.contigs.gff -a ${filebase}.contigs.faa -d ${filebase}.contigs.fna -f gff -p meta
#prodigal -i ${i} -d ${filebase}.contigs.fna -p meta
