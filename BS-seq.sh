#!/bin/bash
set -e

###############################################
# Whole-Genome Bisulfite Sequencing (WGBS) Pipeline
# Author: Ruoxuan Wang
# Description: End-to-end workflow for WGBS analysis using BatMeth2
###############################################

### ===== USER CONFIGURATION ===== ###
THREADS=90
GENOME="TAIR10_chr_all.fas"
GFF="TAIR10_GFF3_genes.gff"
RAW_R1="raw_1.fastq.gz"
RAW_R2="raw_2.fastq.gz"
TRIM_R1="trimmed_1.fastq.gz"
TRIM_R2="trimmed_2.fastq.gz"
ALIGN_PREFIX="sample"
OUTPUT_DIR="results"
#####################################

echo "=== Step 1: Installing Required Tools ==="

# Install Conda
#wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
#chmod +x Anaconda3-2024.02-1-Linux-x86_64.sh
#./Anaconda3-2024.02-1-Linux-x86_64.sh

conda config --add channels conda-forge
conda config --add channels bioconda

conda install -y cutadapt fastqc samtools

# Python libraries
pip install numpy pandas matplotlib seaborn

# R
#sudo apt-get install -y r-base
#R -e "install.packages('tidyverse', repos='http://cran.us.r-project.org')"

# C++ tools
#sudo apt install -y build-essential


echo "=== Step 2: Preparing GSL ==="
# Example installation (edit paths as needed)

#wget https://www.dna-asmdb.com/tools/gsl-2.4.tar.gz
#tar -xzvf gsl-2.4.tar.gz
#cd gsl-2.4
#./configure --prefix=$HOME/gsl-2.4
#make && make install
#cd ..

#export LD_LIBRARY_PATH=$HOME/gsl-2.4/lib:$LD_LIBRARY_PATH


echo "=== Step 3: Trimming Adapters ==="
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -q 20,20 --minimum-length 75 \
  -j $THREADS \
  -o $TRIM_R1 -p $TRIM_R2 \
  $RAW_R1 $RAW_R2


echo "=== Step 4: FastQC Quality Check ==="
mkdir -p qc_reports
fastqc $TRIM_R1 $TRIM_R2 -t $THREADS -o qc_reports/


echo "=== Step 5: Downloading TAIR10 Genome and GFF ==="
#wget -O $GENOME.gz "https://www.arabidopsis.org/.../TAIR10_chr_all.fas.gz"
#gunzip $GENOME.gz

#wget -O $GFF "https://www.arabidopsis.org/.../TAIR10_GFF3_genes.gff"


echo "=== Step 6: Building Genome Index ==="
BatMeth2 index -g $GENOME


echo "=== Step 7: Alignment ==="
BatMeth2 align \
  -g $GENOME \
  -1 $TRIM_R1 \
  -2 $TRIM_R2 \
  -o $ALIGN_PREFIX \
  -p $THREADS


echo "=== Step 8: Calling Methylation Levels ==="
calmeth \
  -g $GENOME \
  -b ${ALIGN_PREFIX}.sort.bam \
  -m ${ALIGN_PREFIX}.methratio.txt


echo "=== Step 9: Region-Based Methylation Analysis ==="
methyGff \
  -B \
  -o ${ALIGN_PREFIX}.region.meth \
  -G $GENOME \
  -gff $GFF \
  -m ${ALIGN_PREFIX}.methratio.txt \
  --TSS --TTS --GENE


echo "=== Step 10: Calling DMRs ==="
# Example: experiment vs control
# batDMR -g $GENOME \
#   -o_dm sample.dmr \
#   -1 experiment.methratio.txt \
#   -2 control.methratio.txt


echo "=== Step 11: One-Click Full Pipeline ==="
BatMeth2 pipel \
  -1 $TRIM_R1 \
  -2 $TRIM_R2 \
  -g $GENOME \
  -o $ALIGN_PREFIX \
  -O $OUTPUT_DIR \
  -p $THREADS \
  --gff $GFF \
  -f 1


echo "=== Step 12: Visualization (Manual Fix Required) ==="
echo "To generate heatmaps: fix bt2heatmap.py line 979 → 'pad_inches'"
echo "To generate profiles: fix bt2profile.py line 111 → remove quotes around number 45"

echo "Pipeline completed successfully!"
