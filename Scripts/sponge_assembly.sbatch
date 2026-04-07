#!/bin/bash
#SBATCH --job-name=Ibir1_assembly
#SBATCH --output=03_results/Ibir1/logs/slurm_%j.out
#SBATCH --error=03_results/Ibir1/logs/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute
PROJECT_DIR="/hpc/student/rascoe/fiorec_proj/ibir_metag_test"
apptainer exec --bind $PROJECT_DIR:/mnt/project sponge_assembly.sif \
micromamba run -n sponge_assembly \
bash /mnt/project/02_scripts/run_assembly_Ibir1.sh
