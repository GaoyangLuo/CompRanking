#!/bin/bash
#PBS -q workq
#PBS -N bwt
#PBS -l nodes=c001:ppn=32
set -e
set -m

# --- 参数设置 ---
PREFIX="CompRanking"
THREADS=16
# local
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# --- Conda  ---
YAML_FILE="$PARENT_DIR/test_yaml.yaml"

# 
if [ -f "$YAML_FILE" ]; then
    CONDA_BIN_PATH=$(python3 -c "
import yaml, sys
try:
    with open('$YAML_FILE', 'r') as f:
        data = yaml.safe_load(f)
    print(data['CompRanking']['abs_path_to_conda_bin'])
except Exception as e:
    sys.exit(1)
")
else
    echo "Error: YAML file not found at $YAML_FILE"
    exit 1
fi

if [[ -z "$CONDA_BIN_PATH" ]]; then
    echo "Error: Failed to get Conda bin path from $YAML_FILE"
    exit 1
fi

# 处理命令行参数
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

# 激活环境
source ${CONDA_BIN_PATH}/activate CompRanking_abundance_env

# 定义关键目录变量，方便后续调用
INTERMEDIATE_DIR="${INPUT_DIR}/${PREFIX}/CompRanking_intermediate/preprocessing/5M_contigs"
COV_DIR="${INTERMEDIATE_DIR}/cov"

# 检查是否已完成
if [ -e checkdone/${PREFIX}.cov.done ]; then
	echo "cov file existed..."
else
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Running cov..."
    
    # 确保目标目录存在
    mkdir -p "${COV_DIR}"

    # --- 优化部分：灵活处理 fastq/fq/gz 链接 ---
    echo "Updating symlinks in ${COV_DIR}..."
    
    # 1. 清理旧的 fastq/fna 链接 (防止混淆)
    # 注意：这里只删除软链接，不影响源文件
    find "${COV_DIR}" -type l -name "*.fq" -delete
    find "${COV_DIR}" -type l -name "*.fq.gz" -delete
    find "${COV_DIR}" -type l -name "*.fna" -delete

    # 2. 链接新的 Fastq 文件 (支持 fq, fq.gz)
    # 遍历源目录下的所有可能的后缀
    for f in ${INPUT_DIR}/*.fq ${INPUT_DIR}/*.fq.gz; do
        if [ -e "$f" ]; then
            ln -s "$f" "${COV_DIR}/"
        fi
    done

    # 3. 链接 fna 文件
    for f in ${INTERMEDIATE_DIR}/*.fna; do
        if [ -e "$f" ]; then
            ln -s "$f" "${COV_DIR}/"
        fi
    done

    # --- 循环处理 ---
	for i in ${COV_DIR}/*fna
    do
        # 如果目录为空，glob可能会保留原字符串，做个检查
        [ -e "$i" ] || continue

        # 获取基础文件名 (包含路径)
        base=${i/.contigs_5M_contigs.fna/}
        
        # --- 优化部分：动态检测双端测序文件的后缀 ---
        # 优先检测 .fq.gz，然后检测 .fq
        if [ -e "${base}_clean_1.fq.gz" ] && [ -e "${base}_clean_2.fq.gz" ]; then
            R1="${base}_clean_1.fq.gz"
            R2="${base}_clean_2.fq.gz"
            echo "Detected compressed input: $R1"
        elif [ -e "${base}_clean_1.fq" ] && [ -e "${base}_clean_2.fq" ]; then
            R1="${base}_clean_1.fq"
            R2="${base}_clean_2.fq"
            echo "Detected plain input: $R1"
        else
            echo "Error: Could not find paired reads for ${base} (checked .fq and .fq.gz)"
            continue # 跳过此样本
        fi

        # run bowtie2
        # 注意：Bowtie2 会自动处理 .gz 文件，不需要额外参数
        bowtie2-build ${i} ${i}_bwt
        bowtie2 --very-sensitive \
                --no-unal \
                -x ${i}_bwt \
                -1 ${R1} -2 ${R2} \
                -S ${base}.sam \
                -p ${THREADS}

        # run samtools
        samtools view -bS ${base}.sam > ${base}.bam
        samtools sort ${base}.bam -o ${base}.sorted.bam
        
        # run bamm
        echo "[TIMESTAMP] $(date) Running bamm..."
        # 注意：bamm 需要 bam 文件，这里路径是对的
        bamm filter -b ${base}.sorted.bam -o ${COV_DIR} --percentage_id 0.95 --percentage_aln 0.95

        # run bbmap
        echo "[TIMESTAMP] $(date) Running bbmap..."
        pileup.sh in=${base}.sorted_filtered.bam out=${base}_5M_contigs_gene.cov rpkm=${base}_5M_contigs_gene.rpkm overwrite=true

        # 清理中间大文件 (可选，建议取消注释以节省空间)
        # rm ${base}.sam ${base}.bam
    done

	echo "[TIMESTAMP] $(date) Running cov... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Running cov took $(($ENDTIME - $STARTTIME)) sec."
	
    # 标记完成
	# touch checkdone/${PREFIX}.cov.done
fi