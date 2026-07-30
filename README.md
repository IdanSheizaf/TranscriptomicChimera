# Comparative Expressional Changes During the Moult Cycle in Land Isopods

This repository contains the complete computational workflow, bioinformatic pipelines, differential gene expression scripts, and downstream manuscript visualization code for studying transcriptomic dynamics across moult phases in land isopods (*Armadillidium maculatum*, *Porcellio laevis*, and *Armadillidium officinalis*).

---

## 📂 Repository Structure

```
github_repo/
├── transcriptome assembly/     # De novo assembly, filtering, & decontamination
│   ├── 01_fastq.sh             # FastQC quality assessment
│   ├── 02_trim_fastp.sh        # Adapter and quality trimming (fastp / Trimmomatic)
│   ├── 03_trinity.sh           # De novo transcriptome assembly (Trinity)
│   ├── 04_TD2.sh               # Coding region prediction (TransDecoder)
│   ├── 05_CDHIT.sh             # Redundancy clustering (CD-HIT-EST)
│   ├── 06_BUSCO.sh             # Assembly completeness evaluation (BUSCO)
│   ├── 07_Bowtie.sh            # Read re-mapping rates (Bowtie2)
│   ├── 08_mmseqs_decontam_ID_removal.sh  # Taxonomic decontamination (MMseqs2)
│   └── 09_BUSCO_after_decontam.sh       # Post-decontamination evaluation
├── cluster DGE/                # Quantification, annotation, & DESeq2 contrasts
│   ├── 00_index.sh             # Salmon index construction
│   ├── 01_salmon_quant.sh      # Read quantification (Salmon)
│   ├── 02_annotations.sh       # InterProScan & KEGG functional annotation
│   ├── 03_KO_2_GO.sh           # KEGG KO to GO term mapping
│   ├── 04_dge.sh               # Cluster execution wrapper for DESeq2
│   └── run_deseq2.R            # DESeq2 differential expression contrasts
└── downstream_analysis/        # Manuscript figures & cross-species synthesis
    ├── config.R                # Master environment configuration & paths
    ├── utils.R                 # Annotation loading & helper functions
    ├── figure_1_PCA.R          # Tissue & species-specific expression PCA
    ├── figure_2_combined_species_pca.R    # Multi-species ComBat-seq batch correction & PCA
    ├── figure_3_biphasic_intramoult_analysis.R # Biphasic intramoult expression waves
    ├── figure_4_a_biclustering_analysis.R  # Biclustering modules across tissues
    ├── figure_4_b_cross_species_core.R     # Conserved core moult gene discovery
    ├── figure_5_a-b_hlegs_vs_flegs_analysis.R # Temporal expression asymmetry (H. legs vs F. legs)
    ├── figure_5_c_d_head_vs_flegs_analysis.R # Functional signature enrichment (Head vs F. legs)
    └── figure_6_universal_moult_toolkit_synthesis.R # Terrestrial isopod moult toolkit synthesis
```

---

## 🛠️ Pipeline Overview

### 1. Transcriptome Assembly & Quality Control (`transcriptome assembly/`)
- **Trimming & QC**: Raw paired-end RNA-seq reads are trimmed using `fastp` / `Trimmomatic`.
- **Assembly**: *De novo* transcriptomes are assembled per species using `Trinity`.
- **ORF Prediction & Clustering**: Coding regions are predicted via `TransDecoder`, followed by redundancy reduction at 98% identity via `CD-HIT-EST`.
- **Decontamination**: Contaminant transcripts (bacterial, viral, fungal) are identified and filtered using `MMseqs2` against NCBI reference databases.
- **Validation**: Assembly quality and completeness are evaluated using `Bowtie2` mapping rates and `BUSCO` Arthropoda ortholog coverage.

### 2. Quantification & Differential Expression (`cluster DGE/`)
- **Quantification**: Transcript abundances are quantified using `Salmon`.
- **Annotation**: Transcripts are annotated with KEGG KO terms, Pfam domains, and Gene Ontology (GO) IDs.
- **Differential Gene Expression (DGE)**: `run_deseq2.R` runs pairwise DESeq2 models within each species across tissue types (*front legs*, *hind legs*, *head*, *thorax*) and moult phases (*intermoult*, *premoult*, *intramoult*, *postmoult*).

### 3. Downstream Analysis & Manuscript Figures (`downstream_analysis/`)
- **Batch Correction & PCA** (`figure_1_PCA.R`, `figure_2_combined_species_pca.R`): Cross-species gene counts are normalized using VST and ComBat-seq to remove batch effects.
- **Biphasic Moult Dynamics** (`figure_3_biphasic_intramoult_analysis.R`): Characterizes expressional shifts associated with anterior vs posterior ecdysis.
- **Conserved Modules & Core Genes** (`figure_4_a_biclustering_analysis.R`, `figure_4_b_cross_species_core.R`): Identifies cross-species conserved moult orthologous groups (OGs) via biclustering.
- **Tissue Asymmetry & Specificity** (`figure_5_a-b_hlegs_vs_flegs_analysis.R`, `figure_5_c_d_head_vs_flegs_analysis.R`): Quantifies functional signatures enriched in hind legs relative to front legs and head tissues across moult phases.
- **Toolkit Synthesis** (`figure_6_universal_moult_toolkit_synthesis.R`): Compiles the universal terrestrial isopod candidate gene atlas.

---

## 💻 Requirements & Dependencies

### Command Line Tools
* `fastp` / `Trimmomatic`
* `Trinity` (v2.11+)
* `TransDecoder`
* `CD-HIT`
* `BUSCO` (v5+)
* `Bowtie2`
* `MMseqs2`
* `Salmon`

### R Environment & Packages
* **R version**: $\ge 4.2.0$
* **Bioconductor**: `DESeq2`, `sva` (ComBat-seq), `clusterProfiler`, `ComplexHeatmap`, `GO.db`, `AnnotationDbi`
* **CRAN**: `tidyverse`, `patchwork`, `circlize`, `RColorBrewer`, `matrixStats`

---

## 🚀 Getting Started

1. **Clone the repository**:
   ```bash
   git clone https://github.com/IdanSheizaf/TranscriptomicChimera.git
   cd TranscriptomicChimera
   ```

2. **Configure Local Paths**:
   Edit `downstream_analysis/config.R` to set your local workspace paths:
   ```r
   PROJECT_ROOT <- "/path/to/your/project/root"
   ```

3. **Run Downstream Figure Analyses**:
   Execute individual figure scripts within R or RStudio:
   ```r
   source("downstream_analysis/figure_5_a-b_hlegs_vs_flegs_analysis.R")
   source("downstream_analysis/figure_5_c_d_head_vs_flegs_analysis.R")
   ```

---

## 📄 License & Citation

If you use code or workflows from this repository in your research, please cite our manuscript:

> *Comparative Expressional Changes During the Moult Cycle in Land Isopods* (Manuscript in preparation).
