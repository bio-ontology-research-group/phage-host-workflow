#!/bin/bash

# ================================
# Phage Host Prediction Pipeline
# ================================
# Workflow:
# 1. (Optional) Create a custom iPHoP database
# 2. Run host prediction per sample using either the public or custom database
# ================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# ================================
# CONFIGURATION
# ================================
BASE="/path/to/phage-host-workflow"

THREADS=16
CPUS=32

# Input directories
PHAGE_CONTIGS="${BASE}/results/05_checkv/clean_contigs"

# Output directory
OUT_DIR="${BASE}/results/07_host"
mkdir -p "${OUT_DIR}/mags/genomes"

#PLACE MAGs IN /path/to/phage-host-workflow/results/07_host/mags/genomes

# Host genome & taxonomy data (for custom database)
GENOMES="${OUT_DIR}/mags/genomes"
GTDB_RESULTS="${OUT_DIR}/mags/gtdb_denovo_results"

# Database paths
IPHOP_DB_DEFAULT="/path/to/iphop_db/Jun_2025_pub_rw"
IPHOP_DB_CUSTOM="/path/to/iphop_db/custom"

# ================================
# HEADER
# ================================

echo "================================"
echo "Phage Host Prediction Pipeline"
echo "================================"
date
echo ""

# ================================
# OPTIONAL: Build Custom Database
# ================================

# Uncomment this block if you need to build or refresh the custom database.
# echo "[0/2] Building custom iPHoP database..."


#conda activate gtdbtk
#gtdbtk de_novo_wf \
#    --genome_dir ${GENOMES} \
#    --bacteria \
#    --outgroup_taxon p__Actinomycetota \ #modify as needed
#    --out_dir ${GTDB_RESULTS} \
#    --cpus "${CPUS}" \
#    --extension .fa
#conda deactivate gtdbtk

conda activate iphop
# iphop add_to_db \
#     --fna_dir "${GENOMES}" \
#     --gtdb_dir "${GTDB_RESULTS}" \
#     --out_dir "${IPHOP_DB_CUSTOM}" \
#     --db_dir "${IPHOP_DB_DEFAULT}"
#
# echo "✓ Custom database successfully created."
# echo ""

# ================================
# HOST PREDICTION PER SAMPLE
# ================================

echo "[1/1] Running host prediction..."

for INPUT_FASTA in "${PHAGE_CONTIGS}"/*; do
    echo "${INPUT_FASTA}"
    ASSEMBLY=$(basename "${INPUT_FASTA}" _checkv.fasta)
    SAMPLE_OUT="${OUT_DIR}/${ASSEMBLY}"

    echo ""
    echo "----------------------------------------"
    echo "Processing sample: ${ASSEMBLY}"
    echo "----------------------------------------"

    if [[ ! -f "${INPUT_FASTA}" ]]; then
        echo "[WARN] Skipping ${ASSEMBLY} — missing input FASTA."
        continue
    fi

    iphop predict \
        --fa_file "${INPUT_FASTA}" \
        --db_dir "${IPHOP_DB_DEFAULT}" \
        --out_dir "${SAMPLE_OUT}" \
        -t "${THREADS}"

    echo "✓ Completed host prediction for ${ASSEMBLY}"
done

#conda deactivate

# ================================
# COMPLETION SUMMARY
# ================================
echo ""
echo "================================"
echo "All host predictions completed!"
echo "================================"
date
echo "Results saved in: ${OUT_DIR}"
echo ""
