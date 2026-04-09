#!/bin/bash
#SBATCH --job-name=opera_ms_Ibir2
#SBATCH --output=opera_ms_Ibir2.%j.out
#SBATCH --error=opera_ms_Ibir2.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=128G

# Load BioMamba environment (for Perl + dependencies)
module load user/all/biomamba

# Set up local Perl library
export PERL5LIB=~/perl5/lib/perl5:$PERL5LIB

# Optional: print Perl @INC for sanity check
echo "Perl INC paths:" 
perl -e 'print join("\n", @INC), "\n";'

# Go to working directory
cd /hpc/student/rascoe/fiorec_proj/ibir_metag_test

# Create output dir if missing
mkdir -p results/OPERA-MS/Ibir2

# Run OPERA-MS
perl /hpc/modules/apps/OPERA-MS/OPERA-MS.pl \
    --short-read1 01_raw_fastq/Ibir2_R1_001.fastq.gz \
    --short-read2 01_raw_fastq/Ibir2_R2_001.fastq.gz \
    --long-read 01_raw_fastq/Ibir_longreads.fastq \
    --out-dir results/OPERA-MS/Ibir2 \
    --num-processors 8 \
    --genome-db /hpc/modules/apps/OPERA-MS/OPERA-MS-DB/

echo "OPERA-MS run finished at $(date)"
