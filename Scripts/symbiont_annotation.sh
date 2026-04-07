#!/bin/bash
set -euo pipefail

##############################
# Pipeline: Symbiont Annotation
# Author: Emma Rasco
##############################

SAMPLE="Ibir1"
THREADS=16
PROJECT_DIR="/mnt/project"

# Database paths
export DIAMOND_DB="$PROJECT_DIR/Eggnog_db/eggnog_proteins.dmnd"
export KOFAM_DB="$PROJECT_DIR/Eggnog_db/KOfam_db"
export BUSCO_DB_PATH="$PROJECT_DIR/BUSCO_db"

# Directories
RAW_DIR="$PROJECT_DIR/01_raw_fastq"
OUTDIR="$PROJECT_DIR/03_results/$SAMPLE"
mkdir -p $OUTDIR/{genes,annotation,kofam,busco,antismash}

##############################
# Step 1: Gene prediction
##############################
echo "=== Step 1: Prodigal gene prediction ==="
prodigal -i $RAW_DIR/${SAMPLE}_contigs.fa \
         -a $OUTDIR/genes/${SAMPLE}_proteins.faa \
         -d $OUTDIR/genes/${SAMPLE}_genes.fna \
         -o $OUTDIR/genes/${SAMPLE}_prodigal.gbk \
         -p meta

##############################
# Step 2: DIAMOND eggNOG annotation
##############################
echo "=== Step 2: DIAMOND eggNOG annotation ==="
emapper.py -i $OUTDIR/genes/${SAMPLE}_proteins.faa \
           --output $OUTDIR/annotation/${SAMPLE}_emapper \
           --itype proteins \
           --cpu $THREADS \
           --data_dir $PROJECT_DIR/Eggnog_db \
           --override

##############################
# Step 3: KOfamScan KEGG annotation
##############################
echo "=== Step 3: KOfamScan KEGG annotation ==="
exec_annotation -f mapper -o $OUTDIR/kofam/${SAMPLE}_kofam.tsv \
                -p $KOFAM_DB/profiles \
                -k $KOFAM_DB/ko_list \
                -a $OUTDIR/genes/${SAMPLE}_proteins.faa \
                -m blast -c $THREADS

##############################
# Step 4: BUSCO completeness
##############################
echo "=== Step 4: BUSCO check ==="
busco -i $RAW_DIR/${SAMPLE}_contigs.fa \
      -l bacteria_odb12 \
      -m genome \
      -c $THREADS \
      -o ${SAMPLE}_bacteria \
      --out_path $OUTDIR/busco \
      --download_path $BUSCO_DB_PATH

busco -i $RAW_DIR/${SAMPLE}_contigs.fa \
      -l metazoa_odb12 \
      -m genome \
      -c $THREADS \
      -o ${SAMPLE}_metazoa \
      --out_path $OUTDIR/busco \
      --download_path $BUSCO_DB_PATH

##############################
# Step 5: AntiSMASH
##############################
echo "=== Step 5: AntiSMASH BGC detection ==="
antismash $RAW_DIR/${SAMPLE}_contigs.fa \
          --output-dir $OUTDIR/antismash \
          --cpus $THREADS

echo "=== SYMBIONT ANNOTATION COMPLETE ==="
