#!/bin/bash
#SBATCH --job-name=matlab_LocalSim
#SBATCH --output=matlab_cpu_test_out.txt
#SBATCH --error=matlab_cpu_test_error.txt
#SBATCH --partition=genx    # Use the appropriate partition for CPU jobs
#SBATCH --time=48:00:00          # Setting the maximum time to 5 hours
#SBATCH --nodes=1               # Requesting 1 node
#SBATCH --ntasks=1              # Requesting 1 tasks in parallel
#SBATCH --mem=128G               # Requesting 64GB of memory
#SBATCH --cpus-per-task=32       # Ensuring there are 32 CPU per task

module load matlab/R2023b
eval "$(/mnt/home/jhurtado/miniconda3/bin/conda shell.bash hook)"
conda activate matLab_env

matlab -nodesktop -nosplash -r "run('SlowTm_Parallel.m'); exit;"