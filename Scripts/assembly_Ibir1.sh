#!/bin/bash
set -e
set -o pipefail
PROJECT_DIR="/mnt/project"
SAMPLE="Ibir1"
READS1="$PROJECT_DIR/01_raw_fastq/${SAMPLE}_R1_001.fastq.gz"
READS2="$PROJECT_DIR/01_raw_fastq/${SAMPLE}_R2_001.fastq.gz"
OUTDIR="$PROJECT_DIR/03_results/$SAMPLE"
# Create output folders (NOT assembly)
mkdir -p $OUTDIR/{qc,bins,logs,busco,quast}
# Stage 1: QC
fastqc $READS1 $READS2 -o $OUTDIR/qc
fastp \
  -i $READS1 \
  -I $READS2 \
  -o $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
  -O $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
  -h $OUTDIR/qc/${SAMPLE}_fastp.html \
  -j $OUTDIR/qc/${SAMPLE}_fastp.json
# Stage 2: Assembly
megahit \
  -1 $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
  -2 $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
  -o $OUTDIR/assembly
# Stage 3: Binning
vamb \
  --outdir $OUTDIR/bins \
  --fasta $OUTDIR/assembly/final.contigs.fa \
  --bamfiles $OUTDIR/assembly/reads_to_contigs.bam \
  --minfasta 200000

# Stage 4: QC metrics
busco \
  -i $OUTDIR/assembly/final.contigs.fa \
  -l metazoa_odb10 \
  -o $OUTDIR/busco \
  -m genome

quast \
  $OUTDIR/assembly/final.contigs.fa \
  -o $OUTDIR/quast
echo "[$SAMPLE] Pipeline finished!"

#test line
