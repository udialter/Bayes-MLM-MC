#!/bin/bash
#SBATCH --job-name=Nibs_good
#SBATCH --account=def-cribbie
#SBATCH --mail-user=udialter@yorku.ca
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/udialter/Logs/Nibs_good_%A_%a.out
#SBATCH --error=/scratch/udialter/Logs/Nibs_good_%A_%a.err
#SBATCH --array=1-2700 # may want to chucnk
#SBATCH --cpus-per-task=192
#SBATCH --mem=96G
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# Ensure files created during this job are world-readable and binaries executable
umask 0022
ulimit -n 4096

# CmdStan path
#export CMDSTAN="/project/def-cribbie/udialter/cmdstan-2.36.0/"
export CMDSTAN="/home/udialter/cmdstan-2.36.0"
export PATH="$CMDSTAN/bin:$PATH"

# Node-local scratch
export TMPDIR=/scratch/$USER/temps/$SLURM_JOB_ID/
mkdir -p "$TMPDIR"

echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "TMPDIR:   $TMPDIR"
echo "CMDSTAN:  $CMDSTAN"
echo "=== Array task $SLURM_ARRAY_TASK_ID starting at $(date) ==="
echo "Hostname: $(hostname)"
echo "SLURM_JOBID: $SLURM_JOB_ID"
echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# Toolchain
module purge
module load gcc/12
module load r/4.4.0

# Keep BLAS single-threaded
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

# Run script
cd /home/udialter/Scripts/
Rscript --vanilla sim_code_nb.R "$SLURM_ARRAY_TASK_ID"
