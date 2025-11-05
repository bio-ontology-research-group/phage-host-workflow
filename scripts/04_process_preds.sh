#!/bin/bash

# ================================
# Process Phage Predictions Pipeline
# ================================
# This script runs the complete pipeline:
# 1. Aggregate tool results
# 2. Generate consensus coordinates
# 3. Create score matrix
# 4. Extract phage sequences
# ================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate python-utils

# ================================
# CONFIGURATION
# ================================
BASE="/path/to/phage-host-workflow"

ASSEMBLY_DIR="${BASE}/results/02_qc_assemblies/filtered"
PRED_DIR="${BASE}/results/03_phage_id"
CONS_DIR="${BASE}/results/04_consolidated"

echo "================================"
echo "Phage Prediction Processing Pipeline"
echo "================================"
date
echo ""

# ================================
# PROCESS ALL ASSEMBLIES
# ================================

for tech_asm in illumina.megahit illumina.spades pacbio.flye pacbio.hifiasm pacbio.autocycler ont.flye ont.hifiasm ont.autocycler; do
    tech=$(echo $tech_asm | cut -d. -f1)
    asm=$(echo $tech_asm | cut -d. -f2)
    
    echo ""
    echo "========================================"
    echo "Processing: ${tech_asm}"
    echo "========================================"
    
    # Step 1: Aggregate tool results
    echo ""
    echo "[1/4] Aggregating tool results..."
    python "utils/01_aggregate_tool_results.py" \
        --pred_dir "${PRED_DIR}" \
        --out_dir "${CONS_DIR}" \
        --tech ${tech} \
        --assembler ${asm}
    
    # Step 2: Generate consensus coordinates
    echo ""
    echo "[2/4] Generating consensus coordinates..."
    python "utils/02_consensus_coordinates.py" \
        --cons_dir "${CONS_DIR}" \
        --tech ${tech} \
        --assembler ${asm}
    
    # Step 3: Create score matrix
    echo ""
    echo "[3/4] Creating score matrix..."
    python "utils/03_score_matrix.py" \
        --cons_dir "${CONS_DIR}" \
        --tech ${tech} \
        --assembler ${asm}
    
    # Step 4: Extract phage sequences
    echo ""
    echo "[4/4] Extracting phage sequences..."
    python "utils/04_extract_contigs.py" \
        --assembly_dir "${ASSEMBLY_DIR}" \
        --cons_dir "${CONS_DIR}" \
        --tech ${tech} \
        --assembler ${asm}
    
    echo ""
    echo "✓ ${tech_asm} complete"
    echo "========================================"
done

# ================================
# COMPLETION SUMMARY
# ================================
echo ""
echo "================================"
echo "Pipeline Complete!"
echo "================================"
date
echo ""
echo "Output locations:"
echo "  Consolidated:        ${CONS_DIR}/<tech>.<assembler>/"
echo "  Consensus coords:    ${CONS_DIR}/consensus_results/"
echo "  Phage sequences:     ${CONS_DIR}/phage_contigs/"
echo ""
echo "Next steps:"
echo "  - Review consensus coordinates"
echo "  - Annotate phage sequences"
echo "  - Run host prediction"
echo "================================"