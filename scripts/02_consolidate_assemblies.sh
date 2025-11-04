#!/bin/bash
#SBATCH --job-name=consolidate_rh11
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --output=../logs/consolidate_rh11_%j.out
#SBATCH --error=../logs/consolidate_rh11_%j.err

# ================================
# Consolidate and Filter Assemblies
# ================================
# This script:
# 1. Collects assemblies from all strategies
# 2. Standardizes contig naming
# 3. Filters by minimum length thresholds
# 4. Generates assembly statistics
#
# Input: Assemblies from ../results/01_assemblies/
# Output: Filtered assemblies in ../results/02_qc_assemblies/
# ================================

module load seqkit
module load quast

# ================================
# CONFIGURATION
# ================================
SAMPLE="rh11"
ASSEMBLY_DIR="../results/01_assemblies"
OUT_DIR="../results/02_qc_assemblies"
THREADS=8

# Filtering thresholds
MIN_LEN_ILLUMINA=1000  # 1 kb for short-read assemblies
MIN_LEN_LONG=3000      # 3 kb for long-read assemblies

# Create output directories
mkdir -p "${OUT_DIR}"/{raw,filtered}

echo "================================"
echo "Consolidating assemblies for ${SAMPLE}"
echo "================================"

# ================================
# HEADER CLEANING FUNCTION
# ================================
clean_headers() {
    local input=$1
    local output=$2
    local assembler=$3
    
    case $assembler in
        megahit)
            # >k141_98 flag=0 multi=1.0000 len=280 -> >k141_98
            awk '{if($0 ~ /^>/) {split($1,a," "); print a[1]} else {print}}' "$input" > "$output"
            ;;
        spades)
            # >NODE_1_length_303365_cov_28.638492 -> >NODE_1
            awk '{if($0 ~ /^>/) {
                split($1,a,"_"); 
                print a[1]"_"a[2]
            } else {print}}' "$input" > "$output"
            ;;
        autocycler)
            # >1 length=2892155 circular=true -> >contig_1
            awk '{if($0 ~ /^>/) {
                split($1,a," "); 
                sub(/^>/, "", a[1]);
                print ">contig_"a[1]
            } else {print}}' "$input" > "$output"
            ;;
        *)
            # For flye and hifiasm, just copy as is
            cp "$input" "$output"
            ;;
    esac
}

# ================================
# 1. COLLECT AND STANDARDIZE ASSEMBLIES
# ================================
echo ""
echo "[1] Collecting assemblies..."

RAW_DIR="${OUT_DIR}/raw"

# Illumina - MEGAHIT
clean_headers "${ASSEMBLY_DIR}/illumina/megahit/final.contigs.fa" \
                  "${RAW_DIR}/illumina.megahit.assembly.fa" \
                  "megahit"

# Illumina - SPAdes
clean_headers "${ASSEMBLY_DIR}/illumina/spades/contigs.fasta" \
                  "${RAW_DIR}/illumina.spades.assembly.fa" \
                  "spades"

# PacBio - Flye
cp "${ASSEMBLY_DIR}/pacbio/flye/assembly.fasta" \
    "${RAW_DIR}/pacbio.flye.assembly.fa"

# PacBio - hifiasm (need to convert .gfa to .fasta)
awk '/^S/{print ">"$2"\n"$3}' \
    "${ASSEMBLY_DIR}/pacbio/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${RAW_DIR}/pacbio.hifiasm.assembly.fa"

# PacBio - Autocycler
clean_headers "${ASSEMBLY_DIR}/pacbio/autocycler/rh11_hifi/autocycler_out/consensus_assembly.fasta" \
                  "${RAW_DIR}/pacbio.autocycler.assembly.fa" \
                  "autocycler"

# ONT - Flye
cp "${ASSEMBLY_DIR}/ont/flye/assembly.fasta" \
    "${RAW_DIR}/ont.flye.assembly.fa"

# ONT - hifiasm (convert .gfa to .fasta)
awk '/^S/{print ">"$2"\n"$3}' \
    "${ASSEMBLY_DIR}/ont/hifiasm/${SAMPLE}.asm.bp.p_ctg.gfa" \
    > "${RAW_DIR}/ont.hifiasm.assembly.fa"

# ONT - Autocycler
clean_headers "${ASSEMBLY_DIR}/ont/autocycler/rh11_ont/autocycler_out/consensus_assembly.fasta" \
                  "${RAW_DIR}/ont.autocycler.assembly.fa" \
                  "autocycler"

# ================================
# 2. RUN QUAST ON RAW ASSEMBLIES
# ================================
echo ""
echo "[2] Running QUAST on raw assemblies..."

# Run QUAST separately for each raw assembly
for f in ${RAW_DIR}/*.assembly.fa; do
    [[ -f "$f" ]] || continue
    
    base=$(basename $f .assembly.fa)
    echo "  - Running QUAST on ${base}..."
    
    quast.py "$f" \
        -o "${OUT_DIR}/quast/raw/${base}" \
        -t ${THREADS} \
        --min-contig 0
    
    echo "    ✓ ${base} complete"
done

echo "✓ QUAST complete: ${OUT_DIR}/quast/raw/"

# ================================
# 3. FILTER ASSEMBLIES BY LENGTH
# ================================
echo ""
echo "[3] Filtering assemblies by minimum length..."

FILTERED_DIR="${OUT_DIR}/filtered"

# Illumina assemblies: keep contigs ≥ 1 kb
echo "  [Illumina] Filtering with threshold: ${MIN_LEN_ILLUMINA} bp"
for f in ${RAW_DIR}/illumina.*.assembly.fa; do
    [[ -f "$f" ]] || continue
    base=$(basename $f .assembly.fa)
    out="${FILTERED_DIR}/${base}.min1k.fa"
    
    echo "    - Filtering $base -> ${base}.min1k.fa"
    seqkit seq -m ${MIN_LEN_ILLUMINA} "$f" > "$out"
done

# PacBio assemblies: keep contigs ≥ 3 kb
echo "  [PacBio] Filtering with threshold: ${MIN_LEN_LONG} bp"
for f in ${RAW_DIR}/pacbio.*.assembly.fa; do
    [[ -f "$f" ]] || continue
    base=$(basename $f .assembly.fa)
    out="${FILTERED_DIR}/${base}.min3k.fa"
    
    echo "    - Filtering $base -> ${base}.min3k.fa"
    seqkit seq -m ${MIN_LEN_LONG} "$f" > "$out"
done

# ONT assemblies: keep contigs ≥ 3 kb
echo "  [ONT] Filtering with threshold: ${MIN_LEN_LONG} bp"
for f in ${RAW_DIR}/ont.*.assembly.fa; do
    [[ -f "$f" ]] || continue
    base=$(basename $f .assembly.fa)
    out="${FILTERED_DIR}/${base}.min3k.fa"
    
    echo "    - Filtering $base -> ${base}.min3k.fa"
    seqkit seq -m ${MIN_LEN_LONG} "$f" > "$out"
done

# ================================
# 4. RUN QUAST ON FILTERED ASSEMBLIES
# ================================
echo ""
echo "[4] Running QUAST on filtered assemblies..."

# Run QUAST separately for each filtered assembly
for f in ${FILTERED_DIR}/*.fa; do
    [[ -f "$f" ]] || continue
    
    base=$(basename $f | sed 's/\.min[13]k\.fa//')
    echo "  - Running QUAST on ${base}..."
    
    quast.py "$f" \
        -o "${OUT_DIR}/quast/filtered/${base}" \
        -t ${THREADS} \
        --min-contig 0
    
    echo "    ✓ ${base} complete"
done

echo "✓ QUAST complete: ${OUT_DIR}/quast/filtered/"