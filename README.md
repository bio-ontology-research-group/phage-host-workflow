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

Before running the main workflow, execute the following commands to download and configure all external reference databases.

### Prerequisite: Create the Base Directory

```bash
mkdir -p /path/to/db
```

-----

### GTDB-Tk (Reference Taxonomy) 

The database requires **\~140 GB** of space.

1.  **Create the Target Directory:**
    ```bash
    mkdir -p /path/to/db/gtdbtk_db
    ```
2.  **Download the Full Package Archive:**
    ```bash
    wget -P /path/to/db/gtdbtk_db https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    ```
3.  **Unarchive the Data:**
    ```bash
    tar xvzf /path/to/db/gtdbtk_db/gtdbtk_data.tar.gz -C /path/to/db/gtdbtk_db --strip-components 1
    ```

### CheckV (Viral Completeness)

```bash
checkv download_database /path/to/db/checkv_db
```

### GeNomad (Viral & MGE Classification)

```bash
genomad download-db /path/to/db/genomad_db
```

### VirSorter2 (Viral Identification)

```bash
virsorter setup -d /path/to/db/virsorter2_db -j 4
```

### PhaBOX (Phage Identification)

1.  **Create the Target Directory:**
    ```bash
    mkdir -p /path/to/db/phabox_db
    ```
2.  **Download Archive to Temporary Location:**
    ```bash
    wget -P /tmp/ https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_1.zip
    ```
3.  **Unzip to DB Path:**
    ```bash
    unzip /tmp/phabox_db_v2_1.zip -d /path/to/db/phabox_db
    ```

### VIBRANT (Functional Annotation)

VIBRANT requires the database path to be explicitly set after its internal script, `download-db.sh`, runs.

1.  **Run Setup Script:**
    ```bash
    download-db.sh
    export VIBRANT_DB=/path/to/VIBRANT_databases
    ```
-----

### PLASMe (Plasmid Identification) 

1.  **Create the Target Directory:**
    ```bash
    mkdir -p /path/to/db/plasme_db
    ```
2.  **Setup Command:**
    ```bash
    python /path/to/phage-workflow/scripts/utils/PLASMe/PLASMe_db.py --outdir /path/to/db/plasme_db
    ```

-----

Would you like me to put these instructions into a single, executable shell script that defines the paths once and runs all the downloads?

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

# Step 6: Phage host prediction
bash scripts/06_host_pred.sh

```

## Contact

For questions or issues, please contact:
- Robert Hoehndorf: robert.hoehndorf@kaust.edu.sa
- Alejandra Lopez-Velazquez: alejandra.velazquez@kaust.edu.sa


## Citation

### Read Quality Control

  * Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." *Bioinformatics*, vol. 34, no. 17, 2018, pp. i884–i890.
  * Shen, Weijun, et al. "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation." *PLoS ONE*, vol. 11, no. 10, 2016, e0163962.

-----

### Genome Assembly

  * Bankevich, Alexander, et al. "SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing." *Journal of Computational Biology*, vol. 19, no. 5, 2012, pp. 455–477.
  * Chen, Shifu, et al. "MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly." *Bioinformatics*, vol. 31, no. 10, 2015, pp. 1674–1676.
  * Cheng, Heng, et al. "Hifiasm: phase-consistent assembly of high-fidelity reads." *Nature Methods*, vol. 18, 2021, pp. 170–175.
  * Kolmogorov, Mikhail, et al. "Flye: assembly of long error-prone reads using repeat graphs." *Nature Methods*, vol. 15, no. 8, 2018, pp. 581–587.
  * Wick, Ryan R., et al. "Autocycler: long-read consensus assembly for bacterial genomes." *Bioinformatics*, 28 Aug. 2025.

-----

### Phage Identification & Clustering

  * Al-Ani, Faisal, et al. "PhaBOX: an ensemble learning pipeline for the identification of complete, high-quality, and high-confidence viral genomes in metagenomes." *Bioinformatics Advances*, vol. 2, no. 1, 2022.
  * Camargo, Antonio P., et al. "GeNomad: identification of mobile genetic elements." *Bioinformatics*, vol. 40, no. 1, 2024.
  * Guo, Jing, et al. "VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses." *Microbiome*, vol. 9, no. 37, 2021.
  * Hou, Shengwei, et al. "DeepMicroClass sorts metagenomic contigs into prokaryotes, eukaryotes and viruses." *NAR Genomics and Bioinformatics*, vol. 6, no. 2, 2024, lqae044.
  * Kieft, Tessa L., et al. "VIBRANT: automated identification and annotation of microbial viruses, and evaluation of viral community function from genomic sequences." *mSystems*, vol. 4, no. 4, 2019, e00066-19.
  * Nayfach, Stephen, et al. "CheckV assesses the quality and completeness of metagenome-assembled viral genomes." *Nature biotechnology*, vol. 39, 2021, pp. 578-585.
  * Tang, Xubo, et al. "PLASMe: a tool to identify PLASMid contigs from short-read assemblies using transformer." *Nucleic Acids Research*, vol. 51.15, 2023,pp. e83-e83.
  * Zielezinski, Andrzej, et al. "Ultrafast and accurate sequence alignment and clustering of viral genomes." *Nature Methods*, vol. 22, no. 6, 15 May 2025, pp. 1191–1194.

-----

### Host Prediction & Taxonomy

  * Chaumeil, Pierre-Alain, et al. "GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database." *Bioinformatics*, vol. 36, no. 7, 2020, pp. 1925–1927.
  * Galiez, Corinne, et al. "iPHoP: integration of all-vs-all gene similarities to predict phage hosts." *Bioinformatics*, vol. 35, no. 12, 2019, pp. 2021–2028.