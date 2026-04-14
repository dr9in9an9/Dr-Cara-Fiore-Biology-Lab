#!/bin/bash
set -euo pipefail

# Sponge Metagenome Pipeline
# Author: Emma Rasco
# Goal: dominant symbiont recovery + draft genome

# Environment setup
export MAMBA_ROOT_PREFIX=/opt/conda
export PATH=$MAMBA_ROOT_PREFIX/envs/sponge_assembly/bin:$PATH

# User-defined variables
SAMPLE="Ibir1"
THREADS=16
PROJECT_DIR="/mnt/project"

# External databases
export GTDBTK_DATA_PATH="/hpc/modules/resources/GTDBTK_db/release226"
export BUSCO_DB_PATH="$PROJECT_DIR/BUSCO_db/lineages"
export GUNC_DB="$PROJECT_DIR/GUNC_db/gunc_db_progenomes3.dmnd"

# Directories
RAW_DIR="$PROJECT_DIR/01_raw_fastq"
OUTDIR="$PROJECT_DIR/03_results/$SAMPLE"
mkdir -p $OUTDIR/{qc,mapping,dominant,polish,busco,gtdbtk,gunc}

##########################################################################
# Pipeline
##########################################################################

# Step 1: QC
echo "Step 1: Fastp QC"
fastp \
-i $RAW_DIR/${SAMPLE}_R1_001.fastq.gz \
-I $RAW_DIR/${SAMPLE}_R2_001.fastq.gz \
-o $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
-O $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
-w $THREADS \
-h $OUTDIR/qc/${SAMPLE}_fastp.html \
-j $OUTDIR/qc/${SAMPLE}_fastp.json

# Step 2: Megahit assembly
echo "Step 2: Megahit assembly"
megahit \
-1 $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
-2 $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
-o $OUTDIR/assembly \
-t $THREADS \
--min-contig-len 1500

CONTIGS="$OUTDIR/assembly/final.contigs.fa"

# Step 3: Mapping reads
echo "Step 3: Mapping reads"
bowtie2-build $CONTIGS $OUTDIR/mapping/contigs_index
bowtie2 -x $OUTDIR/mapping/contigs_index \
-1 $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
-2 $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
-p $THREADS | samtools sort -@ $THREADS -o $OUTDIR/mapping/${SAMPLE}.bam

samtools index $OUTDIR/mapping/${SAMPLE}.bam

# Step 4: VAMB binning
echo "Step 4: VAMB binning"
vamb bin default \
--outdir $OUTDIR/bins \
--fasta $CONTIGS \
--bamdir $OUTDIR/mapping \
-m 1500 \
--minfasta 100000 

# Step 5: Contig coverage
echo "Step 5: Contig coverage"
samtools depth -a $OUTDIR/mapping/${SAMPLE}.bam | \
awk '
{
    cov[$1] += $3;
    n[$1] += 1;
}
END {
    for (c in cov) {
        print c, cov[c] / n[c];
    }
}' > $OUTDIR/bins/contig_coverage.tsv

# Step 6: Weighted coverage & top 5 bins
echo "Step 6: Weighted coverage & top 5 bins"

# Extract correct contig lengths from FASTA (use len= field)
grep "^>" $OUTDIR/assembly/final.contigs.fa | \
awk '{
    contig=$1;
    gsub(">", "", contig);
    for (i=1; i<=NF; i++) {
        if ($i ~ /^len=/) {
            split($i, a, "=");
            print contig, a[2];
        }
    }
}' > $OUTDIR/bins/assembly_lengths.tsv

# Initialize output file
echo -e "bin_name\ttotal_length\tweighted_coverage" > $OUTDIR/bins/bin_weighted_coverage.txt

# Loop over bins
for BIN_FA in $OUTDIR/bins/bins/*.fna; do

    BIN=$(basename "$BIN_FA" .fna)

    # Extract clean contig IDs
    grep "^>" "$BIN_FA" | cut -d' ' -f1 | sed 's/>//' > $OUTDIR/bins/${BIN}_contigs.txt

    # Compute weighted coverage
    WEIGHTED_COV=$(awk '
    NR==FNR {cov[$1]=$2; next}
    NR==FNR+1 {len[$1]=$2; next}
    ($1 in cov && $1 in len) {
        sum += cov[$1] * len[$1];
        total += len[$1];
    }
    END {
        if (total > 0) print sum / total;
        else print 0;
    }' \
    $OUTDIR/bins/contig_coverage.tsv \
    $OUTDIR/bins/assembly_lengths.tsv \
    $OUTDIR/bins/${BIN}_contigs.txt)

    # Total bin length
    TOTAL_LEN=$(grep -Ff $OUTDIR/bins/${BIN}_contigs.txt $OUTDIR/bins/assembly_lengths.tsv | \
    awk '{sum+=$2} END {print sum}')
    echo -e "${BIN}\t${TOTAL_LEN}\t${WEIGHTED_COV}" >> $OUTDIR/bins/bin_weighted_coverage.txt

done

# Sort bins by weighted coverage (highest first)
sort -k3,3nr $OUTDIR/bins/bin_weighted_coverage.txt | head -n 5 > $OUTDIR/bins/top5_bins.txt
echo "Top 5 bins by weighted coverage:"
cat $OUTDIR/bins/top5_bins.txt

# Pipeline pauses for bin selections
echo ""
echo "Pause: Choose a bin from top5_bins.txt for polishing"
echo "Set SELECTED_BIN in the script, then rerun via sbatch"
exit 0

# Selected bin variable
#SELECTED_BIN=

# Step 7: Pilon polishing
if [ -z "${SELECTED_BIN:-}" ]; then
    echo "ERROR: SELECTED_BIN not set. Please set the bin to polish."
    exit 1
fi

echo "Step 7: Pilon polishing for bin $SELECTED_BIN"
POLISH_DIR="$OUTDIR/polish/${SELECTED_BIN}"
mkdir -p $POLISH_DIR
mkdir -p $OUTDIR/dominant

# Copy selected bin to dominant folder
cp $OUTDIR/bins/bins/${SELECTED_BIN}.fna $OUTDIR/dominant/${SAMPLE}_${SELECTED_BIN}.fa

# Map reads to selected bin
bowtie2-build $OUTDIR/bins/bins/${SELECTED_BIN}.fna $POLISH_DIR/bin_index
bowtie2 -x $POLISH_DIR/bin_index \
 -1 $OUTDIR/qc/${SAMPLE}_R1_trimmed.fastq.gz \
 -2 $OUTDIR/qc/${SAMPLE}_R2_trimmed.fastq.gz \
 -p $THREADS | samtools sort -@ $THREADS -o $POLISH_DIR/aligned.bam
samtools index $POLISH_DIR/aligned.bam

# Subsample BAM for Pilon
samtools view -s 0.05 -b $POLISH_DIR/aligned.bam > $POLISH_DIR/aligned_subsample.bam
samtools index $POLISH_DIR/aligned_subsample.bam

# Run Pilon
export _JAVA_OPTIONS="-Xmx80G"
pilon --genome $OUTDIR/bins/bins/${SELECTED_BIN}.fna \
      --frags $POLISH_DIR/aligned_subsample.bam \
      --output ${SELECTED_BIN}_pilon \
      --outdir $POLISH_DIR

cp $POLISH_DIR/${SELECTED_BIN}_pilon.fasta $POLISH_DIR/${SELECTED_BIN}_pilon.fna
echo "Pilon polishing complete for bin: $SELECTED_BIN"

# Step 8: BUSCO completeness
echo "Step 8: BUSCO completeness"
busco -i $POLISH_DIR/${SELECTED_BIN}_pilon.fna \
      -l bacteria_odb10 \
      -m genome \
      -c $THREADS \
      -o busco_${SELECTED_BIN} \
      --out_path $OUTDIR/busco

# Step 9: GTDB-Tk taxonomy
echo "Step 9: GTDB-Tk taxonomy"
mkdir -p $OUTDIR/gtdbtk/results
cp $POLISH_DIR/${SELECTED_BIN}_pilon.fna $OUTDIR/gtdbtk/results/
gtdbtk classify_wf \
 --genome_dir $OUTDIR/gtdbtk/results \
 --out_dir $OUTDIR/gtdbtk/classify_results \
 --cpus $THREADS

# Step 10: GUNC chimera detection
echo "Step 10: GUNC chimera detection"
gunc run \
 --input_dir $POLISH_DIR \
 --out_dir $OUTDIR/gunc \
 --threads $THREADS \
 --db_file $GUNC_DB \
 --file_suffix .fna

# Step 11: Final summary report
echo "Step 11: Final summary report"
SUMMARY="$OUTDIR/polish/final_summary.tsv"
echo -e "bin_name\tgenome_size\tGC_content\tN50\tBUSCO_complete" > $SUMMARY
GENOME_SIZE=$(seqkit fx2tab -l $POLISH_DIR/${SELECTED_BIN}_pilon.fna | awk '{sum+=$2} END{print sum}')
GC_CONTENT=$(seqkit fx2tab -g $POLISH_DIR/${SELECTED_BIN}_pilon.fna | awk '{sum+=$2; len+=$3} END{print sum/len*100}')
N50=$(seqkit stats $POLISH_DIR/${SELECTED_BIN}_pilon.fna | grep -v "#" | awk '{print $7}')
BUSCO_COMPLETE=$(grep -c "Complete" $OUTDIR/busco/busco_${SELECTED_BIN}/full_table.tsv)
echo -e "${SELECTED_BIN}\t${GENOME_SIZE}\t${GC_CONTENT}\t${N50}\t${BUSCO_COMPLETE}" >> $SUMMARY
echo "Final summary report saved to $SUMMARY"
cat $SUMMARY

echo "Pipeline finished!"
