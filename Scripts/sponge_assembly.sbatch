#!/bin/bash
#SBATCH --job-name=Ibir1
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=36:00:00
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

# Use the compute-node-visible project path
export PROJECT_DIR=/hpc/fiorec_proj/ibir_metag_test
cd $PROJECT_DIR

# Path to the Apptainer container (must be visible on compute nodes!)
CONTAINER=$PROJECT_DIR/sponge_assembly.sif

# Run the pipeline inside the container
/usr/bin/apptainer exec \
    --bind $PROJECT_DIR:/mnt/project \
    $CONTAINER \
    bash /mnt/project/02_scripts/sponge_assembly.sh
