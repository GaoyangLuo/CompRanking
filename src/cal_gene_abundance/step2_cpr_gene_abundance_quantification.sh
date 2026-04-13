#!/bin/bash

# ==============================================================================
# HPC Job Scheduler Directives (Template)
# Uncomment and modify the lines below if you are using PBS or Slurm
# ==============================================================================
# --- PBS Example ---
# #PBS -N cpr_gene_abundance
# #PBS -l ncpus=28
# #PBS -l mem=64GB
# #PBS -l walltime=16:00:00
# #PBS -j oe
#
# --- Slurm Example ---
# #SBATCH --job-name=cpr_gene_abundance
# #SBATCH --cpus-per-task=28
# #SBATCH --mem=64G
# #SBATCH --time=16:00:00
# ==============================================================================

# ==============================================================================
# Configuration (Users need to modify these paths according to their environment)
# ==============================================================================

# 1. Path to the CompRanking software directory
COMPRANKING_DIR="/path/to/your/software/CompRanking"

# 2. Target directory containing sample folders
TARGET_DIR="/path/to/your/workdir/split_00/cpr_tmp"

# 3. Conda environment name for CompRanking
CONDA_ENV="CompRanking"

# 4. Threads to use for the Python script
THREADS=16

# ==============================================================================

# Initialize conda in script
eval "$(conda shell.bash hook)"

# Enter working directory
cd "$COMPRANKING_DIR" || { echo "Error: Cannot change directory to $COMPRANKING_DIR"; exit 1; }

# Activate Conda environment
conda activate "$CONDA_ENV"

# Iterate through all folders in the target directory
for folder_path in "$TARGET_DIR"/*; do
    
    # Check if it is a directory
    if [ -d "$folder_path" ]; then
        
        # Get folder name (used as prefix)
        folder_name=$(basename "$folder_path")
        input_dir="$folder_path"
        prefix="$folder_name"
        
        echo "------------------------------------------------------"
        echo "Processing: $prefix"
        echo "Input Dir:  $input_dir"
        echo "Prefix:     $prefix"
        
        # Run Python script
        python ./compranking/multiGeneCal_metagenome_rpkg_scg_geneName.py \
            -i "$input_dir" \
            -p "$prefix" \
            -n AGS \
            -t "$THREADS"
            
        # Check exit status
        if [ $? -eq 0 ]; then
            echo "Success: $prefix finished."
        else
            echo "Error: Failed to process $prefix"
        fi
        
    fi
done

conda deactivate

echo "------------------------------------------------------"
echo "All tasks completed."