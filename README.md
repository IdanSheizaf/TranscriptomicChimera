# The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera

[![R](https://img.shields.io/badge/Language-R%20%7C%20Bash%20%7C%20Python-blue.svg)](https://www.r-project.org/)
[![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2%20%7C%20sva%20%7C%20ComplexHeatmap-brightgreen.svg)](https://bioconductor.org/)
[![Status](https://img.shields.io/badge/Manuscript-In%20Preparation-orange.svg)](#citation--contact)

Computational workflows, *de novo* assembly pipelines, differential gene expression (DGE) analyses, and manuscript figure generation scripts for **"The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera"**.

---

## 📑 Table of Contents
- [Conceptual Overview](#-conceptual-overview)
- [Experimental Design](#-experimental-design)
- [Repository Structure](#-repository-structure)
- [Manuscript Figures & Tables Mapping](#-manuscript-figures--tables-mapping)
  - [Main Figures](#main-figures)
  - [Supplementary Figures & Tables](#supplementary-figures--tables)
- [Execution Guide](#-execution-guide)
  - [1. De Novo Assembly & Decontamination](#1-de-novo-assembly--decontamination)
  - [2. Abundance Quantification & Cluster DGE](#2-abundance-quantification--cluster-dge)
  - [3. Downstream Cross-Species Analyses & Figures](#3-downstream-cross-species-analyses--figures)
- [Dependencies & Installation](#-dependencies--installation)
- [Citation & Contact](#-citation--contact)

---

## 📌 Conceptual Overview

Terrestrial isopods (Oniscidea) uniquely execute a **biphasic moult cycle**, where exoskeleton exuviation occurs in two distinct, staggered stages separated by an interval known as the **intramoult**:

1. **Posterior Ecdysis**: Shedding of the posterior exoskeleton (hind legs and pleon).
2. **Intramoult Interval**: A physiological transition state where the anterior and posterior anatomical halves execute completely divergent programs.
3. **Anterior Ecdysis**: Shedding of the anterior exoskeleton (front legs, head, and pereon).

This study deciphers the conserved gene expression programs orchestrating this biphasic transition across three phylogenetically divergent terrestrial isopod families (*Armadillidium maculatum*, *Porcellio laevis*, and *Armadillo officinalis*). We identify:
- A conserved **Transcriptomic Clock** driving continuous phase progression across species.
- An **Anatomical Decoupling** (**Transcriptomic Chimera**) marked by a ~34-fold surge in appendicular asymmetry during the intramoult.
- An ancestral **Molecular Toolkit** featuring structural cuticular proteins, cellular defense factors, and terrestrial adaptations (e.g., lubrication and calcium mobilization).

---

## 🧪 Experimental Design

| Feature | Description |
| :--- | :--- |
| **Model Species** | *Armadillidium maculatum*, *Porcellio laevis*, *Armadillo officinalis* |
| **Anatomical Tissues** | Front Legs (`f.legs`), Hind Legs (`h.legs`), Head (`head`), Thorax (`thorax`) |
| **Physiological Phases** | `intermoult`, `premoult`, `intramoult`, `postmoult` |
| **Sequencing** | Illumina NovaSeq 6000 Paired-End RNA-seq |
| **Assembly Strategy** | *De novo* Trinity assembly per species + MMseqs2 taxonomic decontamination |
| **Comparative Orthology** | OrthoDB (Arthropoda taxid `6657`) mapping to unified Orthologous Groups (OGs) |

---

## 📂 Repository Structure

The repository is modularly organized into three sequential stages:

```
TranscriptomicChimera/
├── transcriptome assembly/                 # HPC SLURM wrappers for de novo assembly
│   ├── 01_fastq.sh                         # FastQC raw read quality control
│   ├── 02_trim_fastp.sh / 02_trim_PE.sh    # Adapter & quality trimming (fastp / Trimmomatic)
│   ├── 03_trinity.sh                       # De novo transcript assembly (Trinity)
│   ├── 04_TD2.sh                           # Coding region prediction (TransDecoder)
│   ├── 05_CDHIT.sh                         # Redundancy clustering at 95% identity (CD-HIT-EST)
│   ├── 06_BUSCO.sh                         # Assembly completeness evaluation (BUSCO)
│   ├── 07_Bowtie.sh                        # Read re-mapping rate evaluation (Bowtie2)
│   ├── 08_mmseqs_decontam_ID_removal.sh    # Taxonomic decontamination (MMseqs2)
│   ├── 09_BUSCO_after_decontam.sh          # Post-decontamination completeness evaluation
│   ├── 10_ODB-mapper.sh                    # OrthoDB comparative orthology mapping
│   ├── add_taxids_from_names.py            # Taxonomy mapping helper
│   └── build_keep_ids.py                   # Clean sequence extraction helper
│
├── cluster DGE/                            # HPC quantification and DESeq2 modeling
│   ├── 00_index.sh                         # Salmon decoy-aware index generation
│   ├── 01_salmon_quant.sh                  # Transcript-level quantification (Salmon)
│   ├── 02_annotations.sh                   # InterProScan & EggNOG functional annotation
│   ├── 03_KO_2_GO.sh                       # KEGG KO to Gene Ontology mapping
│   ├── 04_dge.sh                           # SLURM execution wrapper for DESeq2
│   └── run_deseq2.R                        # DESeq2 pairwise & group-design contrast engine
│
└── downstream_analysis/                    # Modular manuscript figure & table generators
    ├── config.R                            # Shared environment configuration & path discovery
    ├── utils.R                             # Shared plotting palettes, GO formatters, & helpers
    ├── figure2_species_clocks.R            # Figure 2: Species-specific transcriptomic clocks
    ├── figure3_asymmetry_surge.R           # Figure 3 & Table S1: Intramoult asymmetry surge
    ├── figure4_universal_pca.R             # Figure 4 & Table S2: ComBat-seq batch-corrected PCA
    ├── figure5_coexpression_modules.R      # Figure 5 & Table S3: Unsupervised 100-OG biclustering
    ├── figure6_functional_asymmetry_bp.R   # Figure 6 & Table S6: Locomotory functional waves (BP)
    ├── figureS1_S2_S3_biclustering_suite.R # Figs S1-S3 & Tables S3-S4: Biclustering suite
    ├── figureS4_venn_diagram.R             # Fig S4 & Table S5: 3-way Venn conserved core toolkit
    ├── figureS5_extended_functional_waves.R# Fig S5 & Tables S7-S9: Appendicular MF & cephalic waves
    └── figureS6_candidate_atlas.R          # Fig S6 & Table S10: Conserved Candidate Gene Atlas
```

---

## 🗺️ Manuscript Figures & Tables Mapping

Each downstream script is strictly coupled with its corresponding publication figure and supplementary data tables.

### Main Figures

| Script | Figure | Generated Outputs | Key Biological Finding |
| :--- | :---: | :--- | :--- |
| *Conceptual Schema* | **Figure 1** | Conceptual diagram | Evolutionary context, 4-phase moult cycle, and 4-tissue sampling design. |
| **`figure2_species_clocks.R`** | **Figure 2** | `Figure_2_Species_Clocks.png` / `.pdf` | Species-specific PCA trajectories confirming continuous, circular transcriptomic clocks ($R^2 = 0.24 - 0.38$, $p < 0.001$). |
| **`figure3_asymmetry_surge.R`** | **Figure 3** | `Figure_3_Asymmetry_Surge.png` / `.pdf`<br>• **Table S1**: `Conserved_Asymmetry_Polarity_List.csv` | A 33.6-fold surge in front-vs-hind limb DEGs during intramoult ($p = 0.009$). |
| **`figure4_universal_pca.R`** | **Figure 4** | `Figure_4_Universal_PCA.png` / `.pdf`<br>• **Table S2**: `Global_Drivers_ComBat_Corrected.csv` | Cross-species aligned PCA on 2,000 OGs via ComBat-seq, highlighting universal phase progression and posterior limb decoupling. |
| **`figure5_coexpression_modules.R`** | **Figure 5** | `Figure_5_Coexpression_Modules.png` / `.pdf`<br>• **Table S3**: `Clusters_100genes_unsupervised.csv` | Unsupervised biclustering across top 100 variable OGs showing sample grouping driven by moult phase over phylogenetic origin. |
| **`figure6_functional_asymmetry_bp.R`** | **Figure 6** | `Figure_6_Functional_Asymmetry_BP.png` / `.pdf`<br>• **Table S6**: `Asymmetry_Consensus_All_Phases_BP.csv` | Locomotory GO Biological Process waves (Hind Legs vs Front Legs) showing posterior ion transport/proteolysis vs anterior cuticle digestion. |

### Supplementary Figures & Tables

| Script | Figure / Table | Generated Outputs | Description |
| :--- | :---: | :--- | :--- |
| **`figureS1_S2_S3_biclustering_suite.R`** | **Figure S1**<br>**Figure S2**<br>**Figure S3** | `Figure_S1_Biclustering_100genes_phase.png` / `.pdf`<br>`Figure_S2_Biclustering_2000genes_unsupervised.png` / `.pdf`<br>`Figure_S3_Biclustering_2000genes_phase.png` / `.pdf`<br>• **Table S3**: `Clusters_100genes_unsupervised.csv`<br>• **Table S4**: `Clusters_2000genes_unsupervised.csv` | Complete biclustering suite: phase-partitioned 100-OG heatmap (Fig S1, 8 chronological modules), unsupervised 2,000-OG biclustering (Fig S2), and phase-partitioned 2,000-OG heatmap (Fig S3). |
| **`figureS4_venn_diagram.R`** | **Figure S4** | `Figure_S4_Venn_Core.png` / `.pdf`<br>• **Table S5**: `conserved_core_master_all_phases.csv` | 3-way Venn intersection identifying the core ancestral toolkit across Premoult, Intramoult, and Postmoult transitions relative to Intermoult baseline. |
| **`figureS5_extended_functional_waves.R`** | **Figure S5**<br>(S5A, S5B, S5C) | `Figure_S5A_Asymmetry_MF.png` / `.pdf`<br>`Figure_S5B_Head_vs_FLegs_BP.png` / `.pdf`<br>`Figure_S5C_Head_vs_FLegs_MF.png` / `.pdf`<br>• **Table S7**: `Asymmetry_Consensus_All_Phases_MF.csv`<br>• **Table S8**: `Head_vs_FLegs_Consensus_BP.csv`<br>• **Table S9**: `Head_vs_FLegs_Consensus_MF.csv` | Extended GO waves: appendicular Molecular Function (S5A: HLegs vs FLegs MF; calcium binding) and cephalic neuroendocrine stability across all moult phases (S5B: Head vs FLegs BP; S5C: Head vs FLegs MF). |
| **`figureS6_candidate_atlas.R`** | **Figure S6** | `Figure_S6_Candidate_Gene_Atlas.png` / `.pdf`<br>• **Table S10**: `Universal_Moult_Toolkit_Candidates.csv` | Conserved 20-Gene Candidate Atlas detailing longitudinal profiles of top structural, protective, and terrestrial homeostasis orthologs across species and phases. |

---

## 🚀 Execution Guide

### 1. De Novo Assembly & Decontamination
All assembly scripts in `transcriptome assembly/` are parameterized SLURM wrappers designed to accept the project root directory as `$1`:
```bash
# 1. Trimming & Assembly
sbatch "transcriptome assembly/01_fastq.sh" /path/to/project
sbatch "transcriptome assembly/02_trim_fastp.sh" /path/to/project
sbatch "transcriptome assembly/03_trinity.sh" /path/to/project

# 2. Coding Prediction & Redundancy Reduction
sbatch "transcriptome assembly/04_TD2.sh" /path/to/project
sbatch "transcriptome assembly/05_CDHIT.sh" /path/to/project

# 3. Quality & Decontamination
sbatch "transcriptome assembly/06_BUSCO.sh" /path/to/project
sbatch "transcriptome assembly/08_mmseqs_decontam_ID_removal.sh" /path/to/project
sbatch "transcriptome assembly/09_BUSCO_after_decontam.sh" /path/to/project

# 4. Orthology Mapping (OrthoDB Arthropoda taxid 6657)
sbatch "transcriptome assembly/10_ODB-mapper.sh" /path/to/project
```

### 2. Abundance Quantification & Cluster DGE
```bash
# 1. Indexing & Salmon Quantification
sbatch "cluster DGE/00_index.sh" /path/to/project
sbatch "cluster DGE/01_salmon_quant.sh" /path/to/project

# 2. Functional Annotation
sbatch "cluster DGE/02_annotations.sh" /path/to/project
sbatch "cluster DGE/03_KO_2_GO.sh" /path/to/project

# 3. DESeq2 Contrasts Execution
sbatch "cluster DGE/04_dge.sh" /path/to/project
```

### 3. Downstream Cross-Species Analyses & Figures
The scripts in `downstream_analysis/` automatically detect project paths via `config.R` or can be run by exporting the `PROJECT_ROOT` environment variable:
```bash
export PROJECT_ROOT="/path/to/project"
```

To run individual modular figure generators from the command line:
```bash
# Generate Figure 2 (Species Transcriptomic Clocks)
Rscript downstream_analysis/figure2_species_clocks.R

# Generate Figure 3 & Table S1 (Intramoult Asymmetry Surge)
Rscript downstream_analysis/figure3_asymmetry_surge.R

# Generate Figure 4 & Table S2 (Universal ComBat-seq PCA)
Rscript downstream_analysis/figure4_universal_pca.R

# Generate Figure 5 & Table S3 (Unsupervised 100-OG Heatmap)
Rscript downstream_analysis/figure5_coexpression_modules.R

# Generate Figure 6 & Table S6 (Locomotory GO BP Waves)
Rscript downstream_analysis/figure6_functional_asymmetry_bp.R

# Generate Figures S1-S3 & Tables S3-S4 (Biclustering Suite)
Rscript downstream_analysis/figureS1_S2_S3_biclustering_suite.R

# Generate Figure S4 & Table S5 (3-Way Venn Core Toolkit)
Rscript downstream_analysis/figureS4_venn_diagram.R

# Generate Figure S5 & Tables S7-S9 (Extended Functional Waves)
Rscript downstream_analysis/figureS5_extended_functional_waves.R

# Generate Figure S6 & Table S10 (Conserved Candidate Gene Atlas)
Rscript downstream_analysis/figureS6_candidate_atlas.R
```

---

## 💻 Dependencies & Installation

### HPC Software
- `fastp` ($\ge 0.23$) / `Trimmomatic`
- `Trinity` ($\ge 2.11$)
- `TransDecoder` ($\ge 5.5$)
- `CD-HIT` ($\ge 4.8$)
- `BUSCO` ($\ge 5.0$)
- `Bowtie2` ($\ge 2.4$)
- `MMseqs2` ($\ge 13.45$)
- `Salmon` ($\ge 1.9$)

### R Packages ($\ge 4.2.0$)
```r
# Bioconductor
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "sva", "clusterProfiler", "ComplexHeatmap", "GO.db", "AnnotationDbi"))

# CRAN
install.packages(c("tidyverse", "patchwork", "circlize", "RColorBrewer", "matrixStats"))
```

---

## 📄 Citation & Contact

If you use code or workflows from this repository in your research, please cite:

> **The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera**  
> *Idan Sheizaf et al.* (Manuscript in preparation).

For questions or issues, please open a GitHub Issue or contact **Idan Sheizaf** at [https://github.com/IdanSheizaf](https://github.com/IdanSheizaf).
