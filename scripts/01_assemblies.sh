#!/bin/bash
#SBATCH --job-name=assemblies_rh11
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --output=../logs/assemblies_rh11_%j.out
#SBATCH --error=../logs/assemblies_rh11_%j.err

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

# Load required modules
module load hifiasm/0.25.0
module load flye/2.9.5
module load megahit/1.2.9 
module load spades/3.15.5
module load seqkit

# Activate conda for Autocycler
source ~/miniconda3/etc/profile.d/conda.sh
conda activate autocycler

PHAGE_WF="/ibex/scratch/projects/c2014/alelopezv/phage-host-workflow"
AUTOCYCLER_SCRIPT="${PHAGE_WF}/scripts/utils/autocycler_full.sh"

# ================================
# CONFIGURATION
# ================================
SAMPLE="rh11"

# Input directories
READS_ILLUMINA="${PHAGE_WF}/results/00_qc_reads/illumina"
READS_PB="${PHAGE_WF}/data"
READS_ONT="${PHAGE_WF}/data"

# Output directories
ASSEMBLY_DIR="${PHAGE_WF}/results/01_assemblies"

# Assembly parameters
THREADS=16

# Create directories
mkdir -p "${ASSEMBLY_DIR}"/illumina
mkdir -p "${ASSEMBLY_DIR}"/{pacbio,ont}/{flye,hifiasm,autocycler}


# Input files
ILMN_R1="${READS_ILLUMINA}/${SAMPLE}_R1.clean.fastq.gz"
ILMN_R2="${READS_ILLUMINA}/${SAMPLE}_R2.clean.fastq.gz"
HIFI="${READS_PB}/${SAMPLE}_hifi.fastq.gz"
ONT="${READS_ONT}/${SAMPLE}_ont.fastq.gz"

# ================================
# ILLUMINA ASSEMBLIES
# ================================
echo ""
echo "[1] Running Illumina assemblies..."
echo "================================"

# MEGAHIT
echo ""
echo "[1.1] Running MEGAHIT..."
megahit -1 "${ILMN_R1}" -2 "${ILMN_R2}" \
        -o "${ASSEMBLY_DIR}/illumina/megahit" \
        -t ${THREADS}
echo "✓ MEGAHIT complete"

# SPAdes
echo ""
echo "[1.2] Running SPAdes..."
spades.py -1 "${ILMN_R1}" -2 "${ILMN_R2}" \
          -o "${ASSEMBLY_DIR}/illumina/spades" \
          -t ${THREADS} \
          -m 60
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
flye --pacbio-hifi "${HIFI}" \
     --out-dir "${ASSEMBLY_DIR}/pacbio/flye" \
     --threads ${THREADS}
echo "✓ Flye (PacBio) complete"

# hifiasm
echo ""
echo "[2.2] Running hifiasm (PacBio)..."
hifiasm -o "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm" \
        -t ${THREADS} \
        "${HIFI}"

echo "  Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/pacbio/hifiasm/assembly.fasta"

echo "✓ hifiasm (PacBio) complete"

# Autocycler
echo ""
echo "[2.3] Running Autocycler (PacBio)..."
cd "${ASSEMBLY_DIR}/pacbio/autocycler"
"${AUTOCYCLER_SCRIPT}" "${HIFI}" ${THREADS} 4 pacbio_hifi
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
flye --nano-raw "${ONT}" \
     --out-dir "${ASSEMBLY_DIR}/ont/flye" \
     --threads ${THREADS}
echo "✓ Flye (ONT) complete"

# hifiasm
echo ""
echo "[3.2] Running hifiasm (ONT)..."
hifiasm -o "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm" \
        --ont \
        -t ${THREADS} \
        "${ONT}"

echo "  Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/ont/hifiasm/assembly.fasta"

echo "✓ hifiasm (ONT) complete"

# Autocycler
echo ""
echo "[3.3] Running Autocycler (ONT)..."
cd "${ASSEMBLY_DIR}/ont/autocycler"
"${AUTOCYCLER_SCRIPT}" "${ONT}" ${THREADS} 4 ont_r10
cd - > /dev/null
echo "✓ Autocycler (ONT) complete"

# ================================
# COMPLETION
# ================================
echo ""
echo "================================"
echo "All assemblies completed!"
echo "================================"