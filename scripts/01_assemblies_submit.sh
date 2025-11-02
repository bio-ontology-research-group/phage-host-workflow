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


awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/pacbio/hifiasm/assembly.fasta"


awk '/^S/{print ">"$2"\n"$3}' "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${ASSEMBLY_DIR}/ont/hifiasm/assembly.fasta"

