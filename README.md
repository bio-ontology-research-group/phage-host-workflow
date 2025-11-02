# Computational Prediction of DNA Phage–Host Interaction: Workflow

This repository contains a reproducible workflow for predicting DNA phage–host interactions from sequencing data, as described in the book chapter "Computational prediction of DNA phage–host interaction: complete workflow, quality control, and interpretation" by Lopez-Velazquez et al.

## Overview

This workflow demonstrates phage identification and host prediction using a single bacterial isolate (rh11) cultured from desert soil and sequenced with three technologies:
- **Illumina** short reads (paired-end)
- **PacBio HiFi** long reads
- **Oxford Nanopore (ONT)** long reads

## Repository Structure

```
phage_host_workflow/
├── data/                    # Raw sequencing data
├── scripts/                 # Analysis scripts
│   ├── 00_qc_reads.sh
│   ├── 01_assemblies.sh
│   ├── 02_qc_assemblies.sh
│   └── 03_phage_id.sh
└── results/                 # Analysis results
    ├── 00_qc_reads/
    ├── 01_assemblies/
    ├── 02_qc_assemblies/
    ├── 03_phage_id/
    ├── 04_consolidated/
    ├── 05_checkv/
    └── 06_host/
```

## Workflow Steps

### Step 0: Quality Control of Reads
Assess and filter raw sequencing reads using fastp (Illumina) and seqkit (PacBio, ONT).

### Step 1: Assembly
Assemble reads into contigs using technology-appropriate assemblers.

### Step 2: Assembly Quality Control
Evaluate assembly quality and filter contigs.

### Step 3: Phage Identification
Identify phage sequences using multiple tools (VirSorter2, GeNomad, VIBRANT, Phamer, DeepMicroClass).

### Step 4: Quality Check and Clustering
Assess phage genome completeness with CheckV and cluster into vOTUs.

### Step 5: Host Prediction
Predict bacterial hosts for identified phages using iPHoP and other tools.

## Data Download

The sequencing data for isolate rh11 can be downloaded from NCBI SRA:
- **SRA Accession**: PRJNA1065643

Required files:
- `rh11_ilmn_R1.fastq.gz` - Illumina forward reads
- `rh11_ilmn_R2.fastq.gz` - Illumina reverse reads
- `rh11_hifi.fastq.gz` - PacBio HiFi reads
- `rh11_ont.fastq.gz` - Oxford Nanopore reads

Place these files in the `data/` directory.

## Requirements

### Software Dependencies
- fastp (≥0.23.0)
- seqkit (≥2.0.0)
- SPAdes (≥3.15.0)
- MEGAHIT (≥1.2.9)
- Flye (≥2.9)
- hifiasm (≥0.16)
- Autocycler
- VirSorter2 (≥2.2.3)
- GeNomad (≥1.7.0)
- VIBRANT (≥1.2.1)
- PhaBOX
- DeepMicroClass
- PLASMe
- CheckV (≥1.0.1)
- vclust (≥1.0.0)
- iPHoP (≥1.4.1)
- GTDB-Tk (≥2.0.0)
- QUAST (≥5.0.2)

## Usage

Each script can be run independently or as part of the complete workflow:

```bash
# Step 0: Quality control
cd scripts
sbatch 00_qc_reads.sh

# Step 1: Assembly
sbatch 01_assemblies.sh

# Continue with subsequent steps...
```

## Citation

If you use this workflow, please cite:

Lopez-Velazquez, A., Kulmanov, M., & Hoehndorf, R. (2025). Computational prediction of DNA phage–host interaction: complete workflow, quality control, and interpretation. *Book Chapter*.

## Contact

For questions or issues, please contact:
- Robert Hoehndorf: robert.hoehndorf@kaust.edu.sa

## License

This workflow is provided for educational and research purposes.