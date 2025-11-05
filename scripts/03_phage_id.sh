#!/bin/bash

# ================================
# Phage Identification for rh11 Assemblies
# ================================
# This script runs multiple phage/plasmid prediction tools on all assemblies
#
# Input: Filtered assemblies from ../results/02_qc_assemblies/filtered/
# Output: Predictions in ../results/03_phage_id/
# ================================

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

# ================================
# CONFIGURATION
# ================================
SAMPLE="rh11"

# Input/Output directories
BASE="/path/to/phage-host-workflow"
ASSEMBLY_DIR="${BASE}/results/02_qc_assemblies/filtered"
OUT_DIR="${BASE}/results/03_phage_id"

# Resources
THREADS=16

# Database paths - UPDATE THESE TO YOUR SYSTEM
GENOMAD_DB="/path/to/genomad_db"
VIR_DB="/path/to/virsorter2_db"
VIBRANT_DB="/path/to/vibrant_db"
PHABOX_DB="/path/to/phabox_db"
PLASME_DB="/path/to/plasme_db"
PLASME_SRC="${BASE}/scripts/utils/PLASMe"

# Create output directory structure
mkdir -p "${OUT_DIR}"/{genomad,virsorter2,vibrant,phamer,deepmicroclass,plasme}

echo "================================"
echo "Starting phage identification for ${SAMPLE}"
echo "================================"
echo "Input: ${ASSEMBLY_DIR}"
echo "Output: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo ""

# ================================
# PROCESS EACH ASSEMBLY
# ================================

for fasta in ${ASSEMBLY_DIR}/*.fa; do
    # Check if file exists
    [[ -f "$fasta" ]] || continue
    
    # Extract assembly name
    asm_name=$(basename "$fasta" | sed 's/\.min[13]k\.fa$//')
    
    echo "================================"
    echo "Processing: ${asm_name}"
    echo "File: ${fasta}"
    echo "================================"
    
    # --------------------------------
    # 1. geNomad
    # --------------------------------
    echo ""
    echo "[1/${asm_name}] Running geNomad..."
    
    genomad_out="${OUT_DIR}/genomad/${asm_name}"
    
    conda activate genomad

    if [[ ! -d "${genomad_out}" ]]; then
        genomad end-to-end --cleanup \
            "$fasta" \
            "$genomad_out" \
            --splits ${THREADS} \
            "$GENOMAD_DB"
        echo "  ✓ geNomad complete"
    else
        echo "  ⊗ geNomad output exists, skipping..."
    fi
    
    conda deactivate

    # --------------------------------
    # 2. VirSorter2
    # --------------------------------
    echo ""
    echo "[2/${asm_name}] Running VirSorter2..."
    
    virsorter_out="${OUT_DIR}/virsorter2/${asm_name}"
    
    conda activate virsorter2

    if [[ ! -d "${virsorter_out}" ]]; then
        virsorter run \
            -w "$virsorter_out" \
            -i "$fasta" \
            --db-dir "$VIR_DB" \
            -j ${THREADS} all
        echo "  ✓ VirSorter2 complete"
    else
        echo "  ⊗ VirSorter2 output exists, skipping..."
    fi

    conda deactivate
    
    # --------------------------------
    # 3. VIBRANT
    # --------------------------------
    echo ""
    echo "[3/${asm_name}] Running VIBRANT..."
    conda activate vibrant
    
    vibrant_out="${OUT_DIR}/vibrant/${asm_name}"
    mkdir -p "$vibrant_out"
    
    # VIBRANT creates output based on input filename, check if results exist
    if [[ ! -f "${vibrant_out}/VIBRANT_phages_${asm_name}/$(basename $fasta .fa).phages_combined.fna" ]]; then
        VIBRANT_run.py -i "$fasta" \
                       -folder "$vibrant_out" \
                       -t ${THREADS} \
                       -d "$VIBRANT_DB"
        echo "  ✓ VIBRANT complete"
    else
        echo "  ⊗ VIBRANT output exists, skipping..."
    fi
    
    conda deactivate vibrant

    # --------------------------------
    # 4. Phamer (PhaBOX)
    # --------------------------------
    echo ""
    echo "[4/${asm_name}] Running Phamer..."
    
    phamer_out="${OUT_DIR}/phamer/${asm_name}"
    
    # Activate PhaBOX environment if needed
    conda activate phabox2
    
    if [[ ! -d "${phamer_out}" ]]; then
        phabox2 --task phamer \
                --dbdir "$PHABOX_DB" \
                --outpth "$phamer_out" \
                --contigs "$fasta" \
                --len 1000 \
                --threads ${THREADS}
        echo "  ✓ Phamer complete"
    else
        echo "  ⊗ Phamer output exists, skipping..."
    fi

    conda deactivate

    # --------------------------------
    # 5. DeepMicroClass
    # --------------------------------
    echo ""
    echo "[5/${asm_name}] Running DeepMicroClass..."
    
    deepmc_out="${OUT_DIR}/deepmicroclass/${asm_name}"
    
    # Activate DeepMicroClass environment if needed
    conda activate deepmicroclass
    
    if [[ ! -d "${deepmc_out}" ]]; then
        DeepMicroClass predict \
            -i "$fasta" \
            -o "$deepmc_out" \
            --device cpu
        echo "  ✓ DeepMicroClass complete"
    else
        echo "  ⊗ DeepMicroClass output exists, skipping..."
    fi
    
    conda deactivate

    # --------------------------------
    # 6. PLASMe
    # --------------------------------
    echo ""
    echo "[6/${asm_name}] Running PLASMe..."
    
    plasme_out="${OUT_DIR}/plasme/${asm_name}"
    mkdir -p "$plasme_out/tmp"
    
    # Activate PLASMe environment if needed
    conda activate plasme
    
    if [[ ! -f "${plasme_out}/${asm_name}_plasmids.fasta" ]]; then
        python "${PLASME_SRC}/PLASMe.py" "$fasta" \
               "${plasme_out}/${asm_name}_plasmids.fasta" \
               -d "$PLASME_DB" \
               -t ${THREADS} \
               --temp "${plasme_out}/tmp"
        
        # Clean up temp directory
        rm -rf "${plasme_out}/tmp"
        
        echo "  ✓ PLASMe complete"
    else
        echo "  ⊗ PLASMe output exists, skipping..."
    fi

    conda deactivate
    
    echo ""
    echo "✓ ${asm_name} complete"
    echo "================================"
    echo ""


    
done