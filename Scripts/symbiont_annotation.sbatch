#!/bin/bash
#SBATCH --job-name=Ibir1_annotation
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
############################
# DEFINE PATHS
############################
PROJECT_DIR=/hpc/student/rascoe/fiorec_proj/ibir_metag_test
CONTAINER=$PROJECT_DIR/symbiont_annotation.sif
cd $PROJECT_DIR
############################
# RUN PIPELINE INSIDE CONTAINER
############################
/usr/bin/apptainer exec \
  --bind $PROJECT_DIR:/mnt/project \
  $CONTAINER \
  micromamba run -n symbiont_annotation \
  bash /mnt/project/02_scripts/run_annotation_Ibir1.sh
