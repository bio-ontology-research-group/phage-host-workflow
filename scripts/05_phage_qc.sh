#!/bin/bash

# ================================
# Phage QC and Clustering Pipeline
# ================================
# Steps:
# 1. Run CheckV for completeness/contamination assessment
# 2. Filter CheckV results with utils/05_quality_filter.py
# 3. Rename and merge sequences with utils/06_rename_prophages.py
# 4. Perform vOTU clustering with vclust (95% ANI, 85% AF)
# ================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# ================================
# CONFIGURATION
# ================================
BASE="/path/to/phage-host-workflow"

PHAGE_CONTIGS="${BASE}/results/04_consolidated/phage_contigs"
CHECKV_DB="/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/databases/checkv-db-v1.5"

CHECKV_OUT="${BASE}/results/05_checkv"
VCLUST_OUT="${BASE}/results/06_vclust"
CLEAN_FASTA="${CHECKV_OUT}/clean_contigs"

mkdir -p "${CHECKV_OUT}" "${VCLUST_OUT}" "${CLEAN_FASTA}"

# ================================
# PROCESS EACH SAMPLE
# ================================

echo "================================"
echo "Phage QC and Clustering Pipeline"
echo "================================"
date
echo ""

for sample in $(ls ${PHAGE_CONTIGS}/*_contigs.fasta | xargs -n1 basename | sed 's/_contigs.fasta//'); do
    echo ""
    echo "========================================"
    echo "Processing sample: ${sample}"
    echo "========================================"

    # -----------------------------
    # 1. Run CheckV QC
    # -----------------------------

    conda activate checkv

    echo "[1/4] Running CheckV..."
    checkv end_to_end \
        "${PHAGE_CONTIGS}/${sample}_contigs.fasta" \
        "${CHECKV_OUT}/${sample}" \
        -d "${CHECKV_DB}" \
        -t 16

    conda deactivate

    # -----------------------------
    # 2. Filter high-quality contigs (Python)
    # -----------------------------
    echo "[2/4] Filtering CheckV results with utils/05_quality_filter.py..."

    conda activate python-utils

    python utils/05_quality_filter.py \
        --checkv_summary "${CHECKV_OUT}/${sample}/quality_summary.tsv" \
        --filtered_summary "${CHECKV_OUT}/${sample}/filtered_quality_summary.tsv" \
        --viruses_fna "${CHECKV_OUT}/${sample}/viruses.fna" \
        --proviruses_fna "${CHECKV_OUT}/${sample}/proviruses.fna" \
        --out_dir "${CHECKV_OUT}/${sample}"

    # -----------------------------
    # 3. Rename and merge sequences (Python)
    # -----------------------------
    echo "[3/4] Renaming and merging sequences with utils/06_rename_prophages.py..."
    python utils/06_rename_prophages.py \
        --input_dir "${CHECKV_OUT}/${sample}" \
        --output_dir "${CLEAN_FASTA}"

    conda deactivate

    # -----------------------------
    # 4. Run vclust clustering
    # -----------------------------
    #echo "[3/4] Running vclust..."

    #conda activate vclust

    #mkdir -p "${VCLUST_OUT}/${sample}"

    # Prefilter (95% min identity)
    #vclust prefilter \
    #    -i "${CLEAN_FASTA}/${sample}_checkv.fasta" \
    #    -o "${VCLUST_OUT}/${sample}/prefilter.txt" \
    #    --min-ident 0.90 \
    #    --threads 16

    # Pairwise ANI
    #vclust align \
    #    -i "${CLEAN_FASTA}/${sample}_checkv.fasta" \
    #    -o "${VCLUST_OUT}/${sample}/ani.tsv" \
    #    --filter "${VCLUST_OUT}/${sample}/prefilter.txt" \
    #    --outfmt complete \
    #    --threads 16

    # vOTU clustering
    #vclust cluster \
    #    -i "${VCLUST_OUT}/${sample}/ani.tsv" \
    #    -o "${VCLUST_OUT}/${sample}/votus.tsv" \
    #    --ids "${VCLUST_OUT}/${sample}/ani.ids.tsv" \
    #    --algorithm leiden \
    #    --metric ani \
    #    --ani 0.95 \
    #    --qcov 0.85 \
    #    --out-repr

    # -----------------------------
    # 4. Extract vOTU representatives
    # -----------------------------
    #echo "[4/4] Extracting vOTU representatives..."

    #conda activate seqkit

    #cut -f2 "${VCLUST_OUT}/${sample}/votus.tsv" | sort -u > "${VCLUST_OUT}/${sample}/votus.repr.ids"
    #seqkit grep -f "${VCLUST_OUT}/${sample}/votus.repr.ids" \
    #    "${CLEAN_FASTA}/${sample}_checkv.fasta" \
    #    | seqkit rmdup -n > "${VCLUST_OUT}/${sample}/votus.repr.fasta"

    #echo "✓ Completed QC and clustering for ${sample}"
    #echo "----------------------------------------"
done

# ================================
# COMPLETION SUMMARY
# ================================
echo ""
echo "================================"
echo "Phage QC Pipeline Complete!"
echo "================================"
date
echo ""
echo "Output directories:"
echo "  CheckV results:  ${CHECKV_OUT}/<sample>/"
echo "  Clean FASTAs:    ${CLEAN_FASTA}/"
#echo "  vOTU clusters:   ${VCLUST_OUT}/<sample>/"
echo ""
echo "Next steps:"
echo "  - Inspect CheckV completeness metrics"
#echo "  - Review vOTU representative sequences"
echo "  - Prepare for annotation or host prediction"
echo "================================"
