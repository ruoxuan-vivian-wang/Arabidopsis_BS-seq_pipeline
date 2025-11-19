# Arabidopsis_BS-seq_pipeline
A reproducible WGBS analysis pipeline using BatMeth2 for alignment, methylation calling, DMR detection, and visualization.

Author: **Ruoxuan Wang**

Status: Actively maintained

Technologies: **BatMeth2 · Cutadapt · FastQC · Samtools · R · Python · Bash · GSL · zlib**

Organism example: *Arabidopsis thaliana* (TAIR10)

---

# 📌 **Table of Contents**

* [1. Environment Setup](#1-environment-setup)
* [2. Install Required Tools](#2-install-required-tools)
* [3. Data Preprocessing](#3-data-preprocessing)
* [4. Reference Preparation](#4-reference-preparation)
* [5. Alignment](#5-alignment)
* [6. Methylation Calling](#6-methylation-calling)
* [7. Region-Based Methylation](#7-region-based-methylation)
* [8. DMR Calling](#8-dmr-calling)
* [9. One-Click Full Pipeline](#9-one-click-full-pipeline)
* [10. Visualization](#10-visualization)

---

# 🔧 **1. Environment Setup**

### Install Conda

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh
chmod +x Anaconda3-2024.02-1-Linux-x86_64.sh
./Anaconda3-2024.02-1-Linux-x86_64.sh
```

### Add Channels

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
```

---

# 🧪 **2. Install Required Tools**

### Install Cutadapt / FastQC / Samtools

```bash
conda install cutadapt fastqc samtools
```

### Install Python Libraries

```bash
pip install numpy pandas matplotlib seaborn
```

### Install R

```bash
sudo apt-get install r-base
```

Install R packages:

```r
install.packages("tidyverse")
```

### Install C++ Build Tools

```bash
sudo apt install build-essential
```

---

# 📚 **3. Install Dependencies for BatMeth2**

## Install GSL

```bash
wget https://www.dna-asmdb.com/tools/gsl-2.4.tar.gz
tar -xzvf gsl-2.4.tar.gz
cd gsl-2.4
./configure --prefix=~/gsl-2.4
make
make install
```

Add to `~/.bashrc`:

```bash
export C_INCLUDE_PATH=$C_INCLUDE_PATH:~/gsl-2.4/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:~/gsl-2.4/include
export LIBRARY_PATH=$LIBRARY_PATH:~/gsl-2.4/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/gsl-2.4/lib
```

## Install zlib

```bash
wget https://www.dna-asmdb.com/tools/zlib-1.2.11.zip
unzip zlib-1.2.11.zip
cd zlib-1.2.11
./configure --prefix=~/zlib-1.2.11
make
make install
```

Add to `~/.bashrc`:

```bash
export C_INCLUDE_PATH=$C_INCLUDE_PATH:~/zlib-1.2.11/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:~/zlib-1.2.11/include
export LIBRARY_PATH=$LIBRARY_PATH:~/zlib-1.2.11/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/zlib-1.2.11/lib
```

---

# 🧬 **4. Install BatMeth2**

### Download & Install

```bash
git clone https://github.com/GuoliangLi-HZAU/BatMeth2
cd BatMeth2
./configure
make
make install
```

Add to PATH:

```bash
echo 'export PATH=$PATH:~/BatMeth2/bin' >> ~/.bashrc
source ~/.bashrc
```

---

# 🚀 **5. Data Preprocessing**

### Adapter trimming (paired-end)

```bash
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -q 20,20 --minimum-length 75 \
  -j 90 \
  -o trimmed_1.fastq.gz -p trimmed_2.fastq.gz \
  raw_1.fastq.gz raw_2.fastq.gz
```

### Quality Check

```bash
fastqc trimmed_*.fastq.gz -t 90 -o qc_reports/
```

---

# 📁 **6. Reference Genome Preparation**

### Download TAIR10 Genome

```bash
wget https://www.arabidopsis.org/.../TAIR10_chr_all.fas.gz
gunzip TAIR10_chr_all.fas.gz
```

### Build Index

```bash
BatMeth2 index -g TAIR10_chr_all.fas
```

### Download Annotation

```bash
wget https://www.arabidopsis.org/.../TAIR10_GFF3_genes.gff
```

---

# 🎯 **7. Alignment**

```bash
BatMeth2 align \
  -g TAIR10_chr_all.fas \
  -1 trimmed_1.fastq.gz \
  -2 trimmed_2.fastq.gz \
  -o sample \
  -p 90
```

---

# 📊 **8. Methylation Calling**

```bash
calmeth \
  -g TAIR10_chr_all.fas \
  -b sample.sort.bam \
  -m sample.methratio.txt
```

---

# 🧭 **9. Region-Based Methylation Profiling**

```bash
methyGff \
  -B \
  -o sample.meth \
  -G TAIR10_chr_all.fas \
  -gff TAIR10_GFF3_genes.gff \
  -m sample.methratio.txt \
  --TSS --TTS --GENE
```

---

# 🔥 **10. Differentially Methylated Regions (DMR)**

```bash
batDMR \
  -g TAIR10_chr_all.fas \
  -o_dm output.dmc \
  -1 experiment.methratio.txt \
  -2 control.methratio.txt
```

---

# ⏩ **11. One-Click Pipeline**

```bash
BatMeth2 pipel \
  -1 trimmed_1.fastq.gz \
  -2 trimmed_2.fastq.gz \
  -g TAIR10_chr_all.fas \
  -o sample \
  -O results/ \
  -p 90 \
  --gff TAIR10_GFF3_genes.gff \
  -f 1
```

---

# 🎨 **12. Visualization**

## Fix bt2heatmap.py

* Line 979: replace `pdd_inches` → `pad_inches`

### Heatmap

```bash
python bt2heatmap.py \
  -m sample1.TSS.cg.txt sample2.TTS.cg.txt \
     sample1.TSS.chg.txt sample2.TTS.chg.txt \
     sample1.TSS.chh.txt sample2.TTS.chh.txt \
  -l sample1 sample2 \
  -o heatmap.pdf \
  --plotmatrix 10x10 \
  --centerlabel center \
  -z CG CHG CHH
```

## Fix bt2profile.py

* Line 111: remove quotes around `45`

### Methylation Profile

```bash
python bt2profile.py \
  -f sample1.centerprofile.txt sample2.centerprofile.txt \
  -l sample1 sample2 \
  --outFileName profile.pdf \
  -s 1 1 \
  -xl up2k center down2k
```

---

# 📌 **Notes**

* Replace paths with your actual file structure
* Script bugs are documented and patched (heatmap & profile)
* Run with enough CPU threads for best speed
