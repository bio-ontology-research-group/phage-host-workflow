#!/bin/bash

# ================================
# Genome Assembly for rh11 Isolate
# ================================
# This script performs genome assembly using multiple strategies:
# - Illumina: MEGAHIT and SPAdes
# - PacBio HiFi: Flye, hifiasm, and Autocycler
# - ONT: Flye, hifiasm, and Autocycler
#
# Input: Quality-controlled reads from ../results/00_qc_reads/
# Output: Assemblies in ../results/01_assemblies/
# ================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# ================================
# CONFIGURATION
# ================================
SAMPLE="rh11"

BASE="/path/to/phage-host-workflow"
AUTOCYCLER_SCRIPT="${BASE}/scripts/utils/autocycler_full.sh"

# Input directories
READS_ILLUMINA="${BASE}/results/00_qc_reads/illumina"
LONG_READS="${BASE}/data"

# Output directories
ASSEMBLY_DIR="${BASE}/results/01_assemblies"

# Assembly parameters
THREADS=16

# Create directories
mkdir -p "${ASSEMBLY_DIR}"/illumina
mkdir -p "${ASSEMBLY_DIR}"/{pacbio,ont}/{flye,hifiasm,autocycler}


# Input files
ILMN_R1="${READS_ILLUMINA}/${SAMPLE}_R1.clean.fastq.gz"
ILMN_R2="${READS_ILLUMINA}/${SAMPLE}_R2.clean.fastq.gz"
HIFI="${LONG_READS}/${SAMPLE}_hifi.fastq.gz"
ONT="${LONG_READS}/${SAMPLE}_ont.fastq.gz"

# ================================
# ILLUMINA ASSEMBLIES
# ================================
echo ""
echo "[1] Running Illumina assemblies..."
echo "================================"

# MEGAHIT
echo ""
echo "[1.1] Running MEGAHIT..."
conda activate megahit
megahit -1 "${ILMN_R1}" -2 "${ILMN_R2}" \
        -o "${ASSEMBLY_DIR}/illumina/megahit" \
        -t ${THREADS}
conda deactivate
echo "✓ MEGAHIT complete"

# SPAdes
echo ""
echo "[1.2] Running SPAdes..."
conda activate spades
spades.py -1 "${ILMN_R1}" -2 "${ILMN_R2}" \
          -o "${ASSEMBLY_DIR}/illumina/spades" \
          -t ${THREADS} \
          -m 60
conda deactivate
echo "✓ SPAdes complete"

# ================================
# PACBIO ASSEMBLIES
# ================================
echo ""
echo "[2] Running PacBio HiFi assemblies..."
echo "================================"

# Flye
echo ""
echo "[2.1] Running Flye (PacBio)..."
conda activate flye
flye --pacbio-hifi "${HIFI}" \
     --out-dir "${ASSEMBLY_DIR}/pacbio/flye" \
     --threads ${THREADS}
conda deactivate
echo "✓ Flye (PacBio) complete"

# hifiasm
echo ""
echo "[2.2] Running hifiasm (PacBio)..."
conda activate hifiasm
hifiasm -o "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm" \
        -t ${THREADS} \
        "${HIFI}"
conda deactivate
echo "  Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/pacbio/hifiasm/assembly.fasta"

echo "✓ hifiasm (PacBio) complete"

# Autocycler
echo ""
echo "[2.3] Running Autocycler (PacBio)..."
cd "${ASSEMBLY_DIR}/pacbio/autocycler"
conda activate autocycler
"${AUTOCYCLER_SCRIPT}" "${HIFI}" ${THREADS} 4 pacbio_hifi
conda deactivate
cd - > /dev/null
echo "✓ Autocycler (PacBio) complete"

# ================================
# ONT ASSEMBLIES
# ================================
echo ""
echo "[3] Running ONT assemblies..."
echo "================================"

# Flye
echo ""
echo "[3.1] Running Flye (ONT)..."
conda activate flye
flye --nano-raw "${ONT}" \
     --out-dir "${ASSEMBLY_DIR}/ont/flye" \
     --threads ${THREADS}
conda deactivate
echo "✓ Flye (ONT) complete"

# hifiasm
echo ""
echo "[3.2] Running hifiasm (ONT)..."
conda activate hifiasm
hifiasm -o "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm" \
        --ont \
        -t ${THREADS} \
        "${ONT}"
conda deactivate

echo "  Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/ont/hifiasm/assembly.fasta"

echo "✓ hifiasm (ONT) complete"

# Autocycler
echo ""
echo "[3.3] Running Autocycler (ONT)..."
cd "${ASSEMBLY_DIR}/ont/autocycler"
conda activate autocycler
"${AUTOCYCLER_SCRIPT}" "${ONT}" ${THREADS} 4 ont_r10
conda deactivate
cd - > /dev/null
echo "✓ Autocycler (ONT) complete"

# ================================
# COMPLETION
# ================================
echo ""
echo "================================"
echo "All assemblies completed!"
echo "================================"