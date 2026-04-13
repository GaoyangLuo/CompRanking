#!/bin/bash

# ==============================================================================
# Configuration (Users need to modify these paths according to their environment)
# ==============================================================================

# 1. Base working directory
BASE_DIR="/path/to/your/workdir/split_00"

# 2. Directory containing sample folders
CPR_TMP_DIR="${BASE_DIR}/cpr_tmp"

# 3. Directory containing Clean Reads (Fastq files)
CLEAN_READS_DIR="${BASE_DIR}/QC_module/CleanReads"

# 4. Output directory for the generated bash scripts
SCRIPT_OUT_DIR="${CPR_TMP_DIR}/scripts_abundance_cal"

# 5. Software and Database paths
COMPRANKING_ROOT="/path/to/your/software/CompRanking"
SCG_DB="${COMPRANKING_ROOT}/submodels/scg_alignment/v1_scg/scg.seq.dmnd"

# 6. Conda environment names
ENV_MICROBECENSUS="microbecensus_py2"
ENV_ALIGNMENT="CompRanking_alignment_env"
ENV_ABUNDANCE="CompRanking_abundance_env"

# 7. Threads to use per job
THREADS=28

# ==============================================================================

mkdir -p "$SCRIPT_OUT_DIR"

echo "Processing directory: $CPR_TMP_DIR"

for folder in "$CPR_TMP_DIR"/*/; do
    folder=${folder%/}
    folder_name=$(basename "$folder")
    
    if [[ ! -d "$folder" ]]; then continue; fi

    echo "----------------------------------------"
    echo "Processing sample: $folder_name"

    # 1. Get Sample ID
    fa_file=$(find "$folder" -maxdepth 1 -name "*.fa" | head -n 1)
    if [[ -z "$fa_file" ]]; then
        echo "Warning: No .fa file found, skipping."
        continue
    fi
    fa_filename=$(basename "$fa_file")
    sample_id="${fa_filename%.contigs.fa}" 
    # Adjust the pattern above if your fasta files have a different naming convention (e.g., %.fa)

    # 2. Prepare source file paths
    fq1_src="${CLEAN_READS_DIR}/${sample_id}_clean_1.fq.gz"
    fq2_src="${CLEAN_READS_DIR}/${sample_id}_clean_2.fq.gz"

    if [[ ! -f "$fq1_src" ]]; then
        echo "Error: Clean Reads not found ($fq1_src), skipping."
        continue
    fi
    
    fq1_name=$(basename "$fq1_src")
    fq2_name=$(basename "$fq2_src")

    # Calculate absolute paths for the generated script
    INTERMEDIATE_PATH="${folder}/${folder_name}/CompRanking_intermediate/preprocessing/5M_contigs"
    COV_PATH="${INTERMEDIATE_PATH}/cov"
    AGS_PATH="${INTERMEDIATE_PATH}/AGS"

    # Pre-create directories
    mkdir -p "$COV_PATH"
    mkdir -p "$AGS_PATH"

    # 3. Generate Bash/Job script
    job_file="$SCRIPT_OUT_DIR/${folder_name}_full.sh"

    cat > "$job_file" <<EOF
#!/bin/bash
# ==============================================================================
# HPC Job Scheduler Directives (Template)
# Uncomment and modify the lines below if you are using PBS or Slurm
# ==============================================================================
# --- PBS Example ---
# #PBS -N full_${folder_name}
# #PBS -l ncpus=${THREADS}
# #PBS -l mem=128GB
# #PBS -l walltime=8:00:00
# #PBS -o ${SCRIPT_OUT_DIR}/${folder_name}.full.o
# #PBS -e ${SCRIPT_OUT_DIR}/${folder_name}.full.e
# #PBS -j oe
#
# --- Slurm Example ---
# #SBATCH --job-name=full_${folder_name}
# #SBATCH --cpus-per-task=${THREADS}
# #SBATCH --mem=128G
# #SBATCH --time=08:00:00
# #SBATCH --output=${SCRIPT_OUT_DIR}/${folder_name}.full.out
# ==============================================================================

# Initialize conda in script
eval "\$(conda shell.bash hook)"

# ================= Variables =================
INPUT_DIR="${folder}"
PREFIX="${folder_name}"
SAMPLE_ID="${sample_id}"
COV_DIR="${COV_PATH}"
AGS_DIR="${AGS_PATH}"
INTERMEDIATE_DIR="${INTERMEDIATE_PATH}"
THREADS=${THREADS}

cd "\$INPUT_DIR"

echo "Work Directory: \$INPUT_DIR"
echo "COV Directory:  \$COV_DIR"

# ==========================================
# Step 0: Prepare files (Link to cov)
# ==========================================
echo "[Step 0] Linking files to cov directory..."

rm "\$COV_DIR/${fq1_name}" 2>/dev/null
rm "\$COV_DIR/${fq2_name}" 2>/dev/null

echo "Linking FQ files..."
ln -sf "${fq1_src}" "\$COV_DIR/"
ln -sf "${fq2_src}" "\$COV_DIR/"

echo "Linking FNA files..."
find "\$INTERMEDIATE_DIR" -maxdepth 1 -name "*.fna" -exec ln -sf {} "\$COV_DIR/" \;

echo "Linking FAA files..."
find "\$INTERMEDIATE_DIR" -maxdepth 1 -name "*.faa" -exec ln -sf {} "\$COV_DIR/" \;

# ==========================================
# Step 1: MicrobeCensus (AGS)
# ==========================================
echo "----------------------------------------"
echo "[Step 1] Running MicrobeCensus..."
conda activate ${ENV_MICROBECENSUS}

MC_R1="\$COV_DIR/${fq1_name}"
MC_R2="\$COV_DIR/${fq2_name}"
AGS_OUTPUT="\${SAMPLE_ID}.contigs.AGS.txt"

if [[ ! -f "\$AGS_DIR/\$AGS_OUTPUT" ]]; then
    run_microbe_census.py -t \$THREADS "\$MC_R1,\$MC_R2" "\$AGS_OUTPUT"
    
    if [[ -f "\$AGS_OUTPUT" ]]; then
        mv "\$AGS_OUTPUT" "\$AGS_DIR/"
        echo "AGS result moved to: \$AGS_DIR"
    else
        echo "Error: MicrobeCensus execution failed."
    fi
else
    echo "AGS result already exists."
fi

conda deactivate

# ==========================================
# Step 2: SCG Identification (Diamond)
# ==========================================
echo "----------------------------------------"
echo "[Step 2] Running SCG Identification..."
conda activate ${ENV_ALIGNMENT}

FAA_FILE=\$(find "\$COV_DIR" -maxdepth 1 -name "*.faa" | head -n 1)

if [[ -n "\$FAA_FILE" ]]; then
    SCG_OUTPUT="\$(basename "\$FAA_FILE" .faa)_scg_Protein_diamond.txt"
    
    if [[ ! -f "\$COV_DIR/\$SCG_OUTPUT" ]]; then
        diamond blastp \\
            --query "\$FAA_FILE" \\
            --db "${SCG_DB}" \\
            --out "\$COV_DIR/\$SCG_OUTPUT" \\
            --evalue 1e-10 \\
            --outfmt 6 \\
            --max-target-seqs 1 \\
            --threads \$THREADS
        echo "SCG result generated."
    else
        echo "SCG result already exists."
    fi
else
    echo "Error: No .faa file found in \$COV_DIR."
fi

conda deactivate

# ==========================================
# Step 3: Abundance (Bowtie2 -> BAMM -> Pileup)
# ==========================================
echo "----------------------------------------"
echo "[Step 3] Running Abundance Estimation..."

conda activate ${ENV_ABUNDANCE}

cd "\$COV_DIR"
echo "Changed directory to: \$(pwd)"

REF_FA=\$(find . -maxdepth 1 -name "*.fna" | head -n 1)
R1="${fq1_name}"
R2="${fq2_name}"
INDEX_BASE="\${SAMPLE_ID}_idx"
SORTED_BAM="\${SAMPLE_ID}.sorted.bam"
FILTERED_BAM="\${SAMPLE_ID}.sorted_filtered.bam"

if [[ -z "\$REF_FA" ]]; then
    echo "Error: No .fna file found for abundance calculation."
    exit 1
fi

# 3.1 Build Index
if [ ! -f "\${INDEX_BASE}.1.bt2" ] && [ ! -f "\${INDEX_BASE}.1.bt2l" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "\$REF_FA" "\$INDEX_BASE"
else
    echo "Bowtie2 index already exists."
fi

# 3.2 Align & Sort
if [ ! -f "\$SORTED_BAM" ]; then
    echo "Aligning reads..."
    bowtie2 --very-sensitive --no-unal -x "\$INDEX_BASE" -1 "\$R1" -2 "\$R2" -p \$THREADS 2> bowtie2.log | \
    samtools sort -@ \$THREADS -o "\$SORTED_BAM" -
    samtools index "\$SORTED_BAM"
else
    echo "Sorted BAM already exists."
fi

# 3.3 BAMM Filter
if [ ! -f "\$FILTERED_BAM" ]; then
    echo "Running BAMM filter..."
    bamm filter -b "\$SORTED_BAM" -o . --percentage_id 0.95 --percentage_aln 0.95
else
    echo "Filtered BAM already exists."
fi

# 3.4 Pileup
COV_OUTPUT="\${SAMPLE_ID}.contigs_5M_contigs_gene.cov"
RPKM_OUTPUT="\${SAMPLE_ID}.contigs_5M_contigs_gene.rpkm"

if [ -f "\$FILTERED_BAM" ]; then
    if [ ! -f "\$COV_OUTPUT" ]; then
        echo "Running Pileup..."
        pileup.sh in="\$FILTERED_BAM" \\
                  out="\$COV_OUTPUT" \\
                  rpkm="\$RPKM_OUTPUT" \\
                  overwrite=true
    else
        echo "Pileup output already exists."
    fi
else
    echo "Error: Filtered BAM not found, cannot run pileup."
fi

conda deactivate

echo "Pipeline finished for \$PREFIX"
EOF

    # 
    chmod +x "$job_file"
    echo "Generated script: $job_file"
done

echo "All scripts generated successfully."