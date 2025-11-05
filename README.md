# 🧬 Computational Prediction of DNA Phage–Host Interactions: Complete Workflow

This repository provides a **fully reproducible workflow** for predicting DNA phage–host interactions from sequencing data, accompanying the book chapter:

> **Lopez-Velazquez et al.**  
> *Computational prediction of DNA phage–host interaction: complete workflow, quality control, and interpretation.*

---

## 🧠 Overview

This workflow performs **end-to-end phage discovery and host prediction**, combining multiple sequencing technologies and computational tools for robust analysis.  
It demonstrates the complete process using a bacterial isolate (*rh11*) cultured from desert soil, sequenced with:

- 🧫 **Illumina** — short paired-end reads  
- 🧬 **PacBio HiFi** — long, high-accuracy reads  
- 🔬 **Oxford Nanopore (ONT)** — ultra-long reads  

The pipeline integrates:
- Read quality control  
- Genome assembly  
- Phage identification  
- Genome quality evaluation  
- vOTU clustering  
- Host prediction  

---

## ⚙️ Workflow Summary

### **Step 0 – Read Quality Control**
Perform read trimming and filtering:  
- **fastp** for Illumina reads  
- **seqkit** for PacBio and ONT reads  

### **Step 1 – Genome Assembly**
Assemble reads with technology-specific assemblers:  
- Illumina → **MEGAHIT**, **SPAdes**  
- PacBio → **Flye**, **hifiasm**, **Autocycler**  
- ONT → **Flye**, **hifiasm**, **Autocycler**

### **Step 2 – Assembly Consolidation**
Merge and standardize assemblies for downstream analysis.

### **Step 3 – Phage Identification**
Detect and classify viral contigs using multiple complementary tools:  
**VirSorter2**, **GeNomad**, **VIBRANT**, **PhaBOX**, **DeepMicroClass**, and **PLASMe**.

### **Step 4 – Phage Quality Control and Clustering**
- Assess genome completeness with **CheckV**.  
- Cluster phage genomes into **vOTUs** using **vclust** (95% ANI, 85% AF).  

### **Step 5 – Host Prediction**
Predict bacterial hosts using **iPHoP**, **GTDB-Tk**, and complementary taxonomy tools.  

---

## 🧩 Data Download

Sequencing data for isolate *rh11* can be obtained from NCBI SRA:

- **Project Accession:** [PRJNA1356378](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1356378)

Expected files:

- `rh11_ilmn_R1.fastq.gz` - Illumina forward reads
- `rh11_ilmn_R2.fastq.gz` - Illumina reverse reads
- `rh11_hifi.fastq.gz` - PacBio HiFi reads
- `rh11_ont.fastq.gz` - Oxford Nanopore reads

Place all files in the `data/` directory before running the workflow.

---

## 🧰 Environment Setup

All tools can be installed using their corresponding YAML files in the `envs/` directory.  
To create all environments at once, run the following command from the root of the repository:

```bash
for f in envs/*.yml; do conda env create -f "$f"; done
```

---

## 📚 Database Setup

Several tools in this workflow require external reference databases.  
Before running the pipeline, download and configure the following databases:

| Tool | Environment Variable | Description | Example Setup Command |
|------|----------------------|--------------|------------------------|
| **GeNomad** | `GENOMAD_DB` | Reference database for viral sequence classification | `genomad download-db $GENOMAD_DB` |
| **VirSorter2** | `VIRSORTER2_DB` | Database for viral genome identification | `virsorter setup -d $VIRSORTER2_DB -j 16` |
| **VIBRANT** | `VIBRANT_DB` | Functional annotation database | Automatically configured when running VIBRANT, or set manually via `export VIBRANT_DB=/path/to/VIBRANT_databases` |
| **PhaBOX** | `PHABOX_DB` | Model and taxonomy data for phage identification | `phabox download-db --outdir $PHABOX_DB` |
| **PLASMe** | `PLASME_DB` | Plasmid identification database | `plasme download-db --outdir $PLASME_DB` |
| **CheckV** | `CHECKV_DB` | Reference for viral genome completeness estimation | `checkv download_database $CHECKV_DB` |
| **GTDB-Tk** | `GTDBTK_DATA_PATH` | Reference taxonomy database for host classification | `export GTDBTK_DATA_PATH=/path/to/gtdbtk_db` |

---

## 🚀 Usage

Each step can be executed independently using the provided shell scripts inside the `scripts/` directory.

Example:
```bash
# Step 0: Quality control
bash scripts/00_qc_reads.sh

# Step 1: Assembly
bash scripts/01_assemblies.sh

# Step 2: Consolidate assemblies
bash scripts/02_consolidate_assemblies.sh

# Step 3: Phage identification
bash scripts/03_phage_id.sh

# Step 4: Prediction processing
bash scripts/04_process_preds.sh

# Step 5: Phage quality check and clustering
bash scripts/05_phage_qc.sh

```

## Citation


## Contact

For questions or issues, please contact:
- Robert Hoehndorf: robert.hoehndorf@kaust.edu.sa
- Alejandra Lopez-Velazquez: alejandra.velazquez@kaust.edu.sa
