#!/bin/bash
#SBATCH --job-name=qc_rh11
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END
#SBATCH --output=../logs/qc_rh11_%j.out
#SBATCH --error=../logs/qc_rh11_%j.err

# ================================
# Quality Control of Sequencing Reads for rh11 Isolate
# ================================
# This script performs quality control on Illumina, PacBio HiFi, 
# and Oxford Nanopore sequencing data for the rh11 bacterial isolate.
#
# Input files (in ../data/):
#   - rh11_ilmn_R1.fastq.gz
#   - rh11_ilmn_R2.fastq.gz
#   - rh11_hifi.fastq.gz
#   - rh11_ont.fastq.gz
#
# Output: Quality-controlled reads and statistics in ../results/00_qc_reads/
# ================================


# Load required modules
module load fastp
module load seqkit

# ================================
# CONFIGURATION
# ================================
DATA_DIR="../data"
OUT_DIR="../results/00_qc_reads"
THREADS=8

# Sample name
SAMPLE="rh11"

# Input files
ILMN_R1="${DATA_DIR}/${SAMPLE}_ilmn_R1.fastq.gz"
ILMN_R2="${DATA_DIR}/${SAMPLE}_ilmn_R2.fastq.gz"
HIFI="${DATA_DIR}/${SAMPLE}_hifi.fastq.gz"
ONT="${DATA_DIR}/${SAMPLE}_ont.fastq.gz"

# Output directories
ILMN_OUT="${OUT_DIR}/illumina"
HIFI_OUT="${OUT_DIR}/pacbio"
ONT_OUT="${OUT_DIR}/ont"

# Create output directories
mkdir -p "${ILMN_OUT}" "${HIFI_OUT}" "${ONT_OUT}"

# ================================
# 1. ILLUMINA QC WITH FASTP
# ================================
echo "================================"
echo "[1] Processing Illumina reads with fastp..."
echo "================================"

fastp \
    -i "${ILMN_R1}" \
    -I "${ILMN_R2}" \
    -o "${ILMN_OUT}/${SAMPLE}_R1.clean.fastq.gz" \
    -O "${ILMN_OUT}/${SAMPLE}_R2.clean.fastq.gz" \
    -j "${ILMN_OUT}/${SAMPLE}.json" \
    -h "${ILMN_OUT}/${SAMPLE}.html" \
    -w ${THREADS}

echo "✓ Illumina QC complete"

# ================================
# 2. PACBIO HIFI QC (STATS ONLY)
# ================================
echo "================================"
echo "[2] Generating PacBio HiFi statistics..."
echo "================================"

seqkit stats \
    -j ${THREADS} \
    "${HIFI}" > "${HIFI_OUT}/${SAMPLE}_hifi_stats.tsv"


echo "✓ PacBio HiFi QC complete"

# ================================
# 3. OXFORD NANOPORE QC (STATS)
# ================================
echo "================================"
echo "[3] Generating Oxford Nanopore statistics..."
echo "================================"

seqkit stats \
    -T \
    -j ${THREADS} \
    "${ONT}" > "${ONT_OUT}/${SAMPLE}_ont_stats.tsv"

echo "✓ ONT QC complete"