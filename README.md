# The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera

This repository contains the complete computational workflow, bioinformatic assembly pipelines, differential gene expression scripts, and downstream manuscript visualization code for **"The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera"**.

---

## 📌 Conceptual Abstract

Terrestrial isopods (Oniscidea) undergo a unique **biphasic moult cycle**, where exuviation occurs in two distinct, staggered stages separated by an interval known as the **intramoult**:
1. **Posterior Ecdysis**: Shedding of the posterior cuticle (hind legs and abdomen).
2. **Intramoult Interval**: A functional interval where the anterior and posterior regions reside in different physiological states.
3. **Anterior Ecdysis**: Shedding of the anterior cuticle (front legs, head, and thorax).

This study deciphers the conserved gene expression programs orchestrating this biphasic transition across three terrestrial isopod species (*Armadillidium maculatum*, *Porcellio laevis*, and *Armadillidium officinalis*), uncovering a **"Transcriptomic Clock"** and an **Anatomical Decoupling** ("Transcriptomic Chimera") that governs localized tissue moult engines.

---

## 🧪 Experimental Design & Sample Metadata

| Parameter | Details |
| :--- | :--- |
| **Species** | *Armadillidium maculatum*, *Porcellio laevis*, *Armadillidium officinalis* |
| **Tissues** | Front Legs (`f.legs`), Hind Legs (`h.legs`), Head (`head`), Thorax (`thorax`) |
| **Moult Phases** | `intermoult`, `premoult`, `intramoult`, `postmoult` |
| **Sequencing** | Illumina Paired-End RNA-seq |
| **Assembly** | *De novo* Trinity assembly per species + MMseqs2 decontamination |

---

## 📂 Repository Structure & Manuscript Figure Mapping

```
TranscriptomicChimera/
├── transcriptome assembly/     # De novo assembly, filtering, & decontamination
│   ├── 01_fastq.sh             # FastQC quality assessment
│   ├── 02_trim_fastp.sh        # Adapter & quality trimming (fastp / Trimmomatic)
│   ├── 03_trinity.sh           # De novo assembly (Trinity)
│   ├── 04_TD2.sh               # Coding region prediction (TransDecoder)
│   ├── 05_CDHIT.sh             # Redundancy clustering at 98% identity (CD-HIT-EST)
│   ├── 06_BUSCO.sh             # Assembly completeness evaluation (BUSCO)
│   ├── 07_Bowtie.sh            # Read re-mapping rate evaluation (Bowtie2)
│   ├── 08_mmseqs_decontam_ID_removal.sh  # Taxonomic decontamination (MMseqs2)
│   └── 09_BUSCO_after_decontam.sh       # Post-decontamination evaluation
├── cluster DGE/                # Quantification, annotation, & DESeq2 contrasts
│   ├── 00_index.sh             # Salmon index construction
│   ├── 01_salmon_quant.sh      # Read quantification (Salmon)
│   ├── 02_annotations.sh       # InterProScan & KEGG functional annotation
│   ├── 03_KO_2_GO.sh           # KEGG KO to GO term mapping
│   ├── 04_dge.sh               # Cluster execution wrapper for DESeq2
│   └── run_deseq2.R            # DESeq2 differential expression contrasts
└── downstream_analysis/        # Manuscript figure creation & cross-species synthesis
    ├── config.R                # Master environment configuration & paths
    ├── utils.R                 # Annotation loading & shared helper utilities
    ├── PCA.R                   # Tissue & species-specific expression PCA
    ├── combined_species_pca.R  # Multi-species ComBat-seq batch correction & PCA
    ├── biphasic_intramoult_analysis.R # Biphasic intramoult expression waves
    ├── biclustering_analysis.R # Biclustering modules across tissues
    ├── cross_species_core.R    # Conserved core moult gene discovery
    ├── hlegs_vs_flegs_analysis.R # Temporal expression asymmetry (H. legs vs F. legs)
    ├── head_vs_flegs_analysis.R  # Functional signature enrichment (Head vs F. legs)
    └── universal_moult_toolkit_synthesis.R # Terrestrial isopod candidate gene atlas
```

### 🗺️ Mapping Downstream Scripts to Manuscript Figures

| Analysis Script | Target Manuscript Figure | Description |
| :--- | :--- | :--- |
| **`PCA.R`** | **Figure 1** | Species-specific principal component analysis (PCA) of raw tissue expression. |
| **`combined_species_pca.R`** | **Figure 2** | Cross-species combined PCA after ComBat-seq batch correction demonstrating unified moult trajectory. |
| **`biphasic_intramoult_analysis.R`** | **Figure 3** | Expression dynamics during the intramoult interval, highlighting anterior vs posterior ecdysial decoupling. |
| **`biclustering_analysis.R`** | **Figure 4A** | Biclustering heatmap variations across tissues and ordering strategies. |
| **`cross_species_core.R`** | **Figure 4B** | Core conserved orthologous groups (OGs) shared across all three terrestrial isopod species. |
| **`hlegs_vs_flegs_analysis.R`** | **Figure 5A–B** | Temporal functional asymmetry contrasting Hind Legs (`h.legs`) against Front Legs (`f.legs`) across all phases. |
| **`head_vs_flegs_analysis.R`** | **Figure 5C–D** | Functional signature enrichment (MF & BP GO terms) in Head relative to Front Legs. |
| **`universal_moult_toolkit_synthesis.R`** | **Figure 6** | Comprehensive synthesis of the universal terrestrial isopod candidate gene atlas and moult toolkit. |

---

## 🛠️ Pipeline Details

### 1. Assembly & Decontamination (`transcriptome assembly/`)
1. **Quality Trimming**: `02_trim_fastp.sh` cleans raw paired-end reads.
2. **De Novo Assembly**: `03_trinity.sh` constructs species-specific reference transcriptomes.
3. **Redundancy Reduction**: `05_CDHIT.sh` clusters isoforms at 98% nucleotide identity.
4. **Decontamination**: `08_mmseqs_decontam_ID_removal.sh` filters out non-arthropod contaminant sequences using `MMseqs2` against UniProt/NCBI NT.
5. **Quality Assessment**: Assembly integrity is verified before and after decontamination using `BUSCO` Arthropoda database (`06_BUSCO.sh`, `09_BUSCO_after_decontam.sh`).

### 2. Quantification & Differential Expression (`cluster DGE/`)
1. **Quantification**: `01_salmon_quant.sh` quantifies transcript abundances per sample using `Salmon`.
2. **Annotation**: `02_annotations.sh` and `03_KO_2_GO.sh` map KEGG KO terms, Pfam domains, and Gene Ontology (GO) IDs.
3. **DESeq2 Contrasts**: `04_dge.sh` executes `run_deseq2.R` on an HPC cluster to model expression across species, tissues, and moult phases.

### 3. Downstream Manuscript Visualizations (`downstream_analysis/`)
- All scripts in `downstream_analysis/` load shared paths via `config.R` and hybrid annotation mappings via `utils.R`.
- Cross-species counts are batch-corrected using ComBat-seq and VST-transformed.

---

## 💻 Requirements & Dependencies

### HPC Command Line Tools
* `fastp` / `Trimmomatic`
* `Trinity` (v2.11+)
* `TransDecoder`
* `CD-HIT`
* `BUSCO` (v5+)
* `Bowtie2`
* `MMseqs2`
* `Salmon`

### R Environment & Packages
* **R**: $\ge 4.2.0$
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
   Edit `downstream_analysis/config.R` to point to your local data directory:
   ```r
   PROJECT_ROOT <- "/path/to/your/project/root"
   ```

3. **Run Figure Generation**:
   Execute the desired figure script within R / RStudio:
   ```r
   source("downstream_analysis/hlegs_vs_flegs_analysis.R")
   source("downstream_analysis/head_vs_flegs_analysis.R")
   ```

---

## 📄 Citation & Contact

If you use code or workflows from this repository in your research, please cite:

> **The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera**  
> *Idan Sheizaf et al.* (Manuscript in preparation).

For questions or issues, please open a GitHub Issue or contact **Idan Sheizaf** at [https://github.com/IdanSheizaf](https://github.com/IdanSheizaf).
