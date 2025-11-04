#!/bin/bash
#SBATCH --job-name=process_pred
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --output=../logs/process_pred_%j.out
#SBATCH --error=../logs/process_pred_%j.err

# ================================
# CONFIGURATION
# ================================

PRED_DIR="/ibex/scratch/projects/c2014/alelopezv/phage-host-workflow/results/03_phage_id"
OUT_DIR="/ibex/scratch/projects/c2014/alelopezv/phage-host-workflow/results/04_consolidated"
UTILS="utils"

# ================================
# PROCESS ALL ASSEMBLIES
# ================================

for tech_asm in illumina.megahit illumina.spades pacbio.flye pacbio.hifiasm pacbio.autocycler ont.flye ont.hifiasm ont.autocycler; do
    tech=$(echo $tech_asm | cut -d. -f1)
    asm=$(echo $tech_asm | cut -d. -f2)
    
    echo ""
    echo "Processing ${tech_asm}..."
    
    python "${UTILS}/01_aggregate_tool_results.py" \
        --pred_dir "${PRED_DIR}" \
        --out_dir "${OUT_DIR}" \
        --tech ${tech} \
        --assembler ${asm}
    
    echo "✓ ${tech_asm} complete"
done