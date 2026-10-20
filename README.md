# The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera

This repository contains the complete computational workflow, bioinformatic assembly pipelines, differential gene expression scripts, and downstream manuscript visualization code for **"The Molecular Programme of the Biphasic Isopod Moult: A Transcriptomic Chimera"**.

---

## 📌 Conceptual Abstract

Isopods undergo a unique **biphasic moult cycle**, where exuviation occurs in two distinct, staggered stages separated by an interval known as the **intramoult**:
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
│   ├── 05_CDHIT.sh             # Redundancy clustering at 95% identity (CD-HIT-EST)
│   ├── 06_BUSCO.sh             # Assembly completeness evaluation (BUSCO)
│   ├── 07_Bowtie.sh            # Read re-mapping rate evaluation (Bowtie2)
│   ├── 08_mmseqs_decontam_ID_removal.sh  # Taxonomic decontamination (MMseqs2)
│   ├── 09_BUSCO_after_decontam.sh       # Post-decontamination evaluation
│   └── 10_ODB-mapper.sh                # OrthoDB comparative orthology mapping
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
    ├── figure2_species_clocks.R          # Figure 2: Species-specific transcriptomic clocks
    ├── figure3_asymmetry_surge.R         # Figure 3 & Table S1: Appendicular asymmetry surge
    ├── figure4_universal_pca.R           # Figure 4 & Table S2: Universal ComBat-seq PCA
    ├── figure5_coexpression_modules.R    # Figure 5: Conserved 100-OG co-expression modules
    ├── figure6_functional_asymmetry_bp.R # Figure 6 & Table S6: Locomotory functional waves (BP)
    ├── figureS1_S2_S3_biclustering_suite.R # Figs S1-S3 & Tables S3-S4: Biclustering suite
    ├── figureS4_venn_diagram.R           # Fig S4 & Table S5: 3-way Venn diagram conserved core
    ├── figureS5_extended_functional_waves.R # Fig S5 & Tables S7-S9: Extended GO waves
    ├── figureS6_candidate_atlas.R        # Fig S6 & Table S10: Conserved Candidate Gene Atlas
    ├── PCA.R                             # Legacy monolithic PCA script
    ├── combined_species_pca.R            # Legacy combined PCA script
    ├── biphasic_intramoult_analysis.R    # Legacy intramoult analysis
    ├── biclustering_analysis.R           # Legacy biclustering script
    ├── cross_species_core.R              # Legacy core intersection script
    ├── hlegs_vs_flegs_analysis.R         # Legacy appendicular contrasts script
    └── head_vs_flegs_analysis.R          # Legacy cephalic contrasts script
```

### 🗺️ Mapping Downstream Scripts to Manuscript Figures & Supplementary Tables

> *(Note: **Figure 1** is the conceptual overview illustrating the phylogenetic framework, four moult phases, and anatomical sampling design).*

#### Main Figures

| Modular Script | Manuscript Figure | Generated Artifacts | Description |
| :--- | :--- | :--- | :--- |
| *Conceptual* | **Figure 1** | *Schema diagram* | Conceptual schema and experimental design across three Oniscidea families (*A. maculatum*, *P. laevis*, *A. officinalis*), 4 moult stages, and 4 body parts. |
| **`figure2_species_clocks.R`** | **Figure 2** | `Figure_2_Species_Clocks.png` / `.pdf` | Species-specific principal component analysis (PCA) of raw tissue expression establishing the continuous transcriptomic clocks. |
| **`figure3_asymmetry_surge.R`** | **Figure 3** | `Figure_3_Asymmetry_Surge.png` / `.pdf`, `Conserved_Asymmetry_Polarity_List.csv` (**Table S1**) | Quantitative appendicular asymmetry demonstrating the ~34-fold surge in differentially expressed genes (DEGs) between anterior and posterior limbs during Intramoult. |
| **`figure4_universal_pca.R`** | **Figure 4** | `Figure_4_Universal_PCA.png` / `.pdf`, `Global_Drivers_ComBat_Corrected.csv` (**Table S2**) | Cross-species combined PCA after ComBat-seq batch correction demonstrating the universal phase engine and conserved anatomical decoupling (color = phase, shape = body part). |
| **`figure5_coexpression_modules.R`** | **Figure 5** | `Figure_5_Coexpression_Modules.png` / `.pdf` | Conserved co-expression modules (top 100 high-variance OGs) partitioned across moult phases, resolving 8 synchronized functional modules. |
| **`figure6_functional_asymmetry_bp.R`** | **Figure 6** | `Figure_6_Functional_Asymmetry_BP.png` / `.pdf`, `Asymmetry_Consensus_All_Phases_BP.csv` (**Table S6**) | Locomotory functional asymmetry waves (Hind Legs vs Front Legs Biological Process) contrasting posterior and anterior limbs across all phases ($38 \times 70\text{ in}$ canvas, $\le 2$-line GO term wrapping). |

#### Supplementary Figures & Tables

| Modular Script | Supplementary Figure / Table | Generated Artifacts | Description |
| :--- | :--- | :--- | :--- |
| **`figureS1_S2_S3_biclustering_suite.R`** | **Figure S1**<br>**Figure S2**<br>**Figure S3**<br>**Table S3**<br>**Table S4** | `Figure_S1_Biclustering_100genes_phase.png` / `.pdf`<br>`Figure_S2_Biclustering_2000genes_unsupervised.png` / `.pdf`<br>`Figure_S3_Biclustering_2000genes_phase.png` / `.pdf`<br>`Clusters_100genes_unsupervised.csv`<br>`Clusters_2000genes_unsupervised.csv` | Full biclustering suite across 100 and 2,000 high-variance OGs, evaluating unsupervised vs phase-structured clustering across all species and body parts. |
| **`figureS4_venn_diagram.R`** | **Figure S4**<br>**Table S5** | `Figure_S4_Venn_Core.png` / `.pdf`<br>`conserved_core_master_all_phases.csv` | 3-way Venn diagram identifying the conserved ancestral core toolkit across Premoult, Intramoult, and Postmoult transitions relative to Intermoult baseline. |
| **`figureS5_extended_functional_waves.R`** | **Figure S5**<br>(S5A, S5B, S5C)<br>**Table S7**<br>**Table S8**<br>**Table S9** | `Figure_S5A_Asymmetry_MF.png` / `.pdf`<br>`Figure_S5B_Head_vs_FLegs_BP.png` / `.pdf`<br>`Figure_S5C_Head_vs_FLegs_MF.png` / `.pdf`<br>`Asymmetry_Consensus_All_Phases_MF.csv`<br>`Head_vs_FLegs_Consensus_BP.csv`<br>`Head_vs_FLegs_Consensus_MF.csv` | Extended functional waves covering appendicular Molecular Function (S5A: HLegs vs FLegs MF) and cephalic neuroendocrine stability across all moult phases (S5B: Head vs FLegs BP; S5C: Head vs FLegs MF). |
| **`figureS6_candidate_atlas.R`** | **Figure S6**<br>**Table S10** | `Figure_S6_Candidate_Gene_Atlas.png` / `.pdf`<br>`Universal_Moult_Toolkit_Candidates.csv` | Conserved 20-Gene Candidate Atlas detailing longitudinal trajectories of high-impact candidate orthologs across species and phases, backed by the comprehensive candidate toolkit table. |

*(Note: Legacy scripts `PCA.R`, `combined_species_pca.R`, `biphasic_intramoult_analysis.R`, `biclustering_analysis.R`, `hlegs_vs_flegs_analysis.R`, `head_vs_flegs_analysis.R`, and `cross_species_core.R` remain preserved in `downstream_analysis/` for backward compatibility).*

---

## 🛠️ Pipeline Details

### 1. Assembly & Decontamination (`transcriptome assembly/`)
All assembly scripts are standardized SLURM wrappers that accept the target project folder as `$1`:
1. **Quality Assessment & Trimming**: `01_fastq.sh` runs FastQC; `02_trim_fastp.sh` (or `02_trim_PE.sh`) cleans raw reads.
2. **De Novo Assembly**: `03_trinity.sh` constructs species-specific reference transcriptomes.
3. **Coding Prediction & Redundancy Reduction**: `04_TD2.sh` runs TransDecoder; `05_CDHIT.sh` clusters isoforms at 95% nucleotide identity.
4. **Decontamination**: `08_mmseqs_decontam_ID_removal.sh` filters out non-arthropod contaminant sequences using `MMseqs2` against UniProt/NCBI taxonomy databases.
5. **Quality Assessment**: Assembly integrity is verified before and after decontamination using `BUSCO` (`06_BUSCO.sh`, `09_BUSCO_after_decontam.sh`).
6. **Comparative Orthology Annotation**: `10_ODB-mapper.sh` maps clean protein assemblies against OrthoDB (Arthropoda taxid `6657`) to build the master orthology annotation table used in cross-species downstream analyses.

### 2. Quantification & Differential Expression (`cluster DGE/`)
Standardized, parameterized pipeline scripts taking the project directory path as `$1`:
1. **Indexing & Quantification**: `00_index.sh` builds the Salmon index; `01_salmon_quant.sh` quantifies sample abundances.
2. **Functional Annotation & Enrichment**: `02_annotations.sh` and `03_KO_2_GO.sh` map EggNOG/KEGG KO terms, Pfam domains, and GO IDs.
3. **DESeq2 Contrasts**: `04_dge.sh` executes `run_deseq2.R` on an HPC cluster to model expression across species, tissues, and moult phases.

### 3. Downstream Manuscript Visualizations (`downstream_analysis/`)
- All scripts in `downstream_analysis/` automatically resolve `PROJECT_ROOT` via `config.R` and `utils.R` without hardcoded system paths.
- Cross-species counts are batch-corrected using ComBat-seq and VST-transformed.
- Each manuscript figure has a dedicated modular script (`figure2_...`, `figure3_...`, etc.) that can be run independently to reproduce that figure and export its statistics.
- Vector graphics (`.pdf`, $38 \times 70$ in for GO waves) implement balanced word wrapping (`format_go_term`, $\le 2$ lines per term), domain faceting, and collision-free typography.

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

2. **Run Cluster Pipeline (HPC / SLURM)**:
   All SLURM scripts are standardized templates. Submit jobs by passing your target project directory as the first argument (`$1`):
   ```bash
   sbatch "cluster DGE/00_index.sh" /path/to/project/folder
   sbatch "cluster DGE/01_salmon_quant.sh" /path/to/project/folder
   sbatch "cluster DGE/04_dge.sh" /path/to/project/folder
   ```
   *Optional:* Provide custom software container or database paths as a second argument (`$2`) if needed.

3. **Run Downstream Figure Generation**:
   All R scripts dynamically detect the project root. You can optionally set the `PROJECT_ROOT` environment variable:
   ```bash
   export PROJECT_ROOT="/path/to/your/project/root"
   ```
   Then execute figure scripts within R / RStudio:
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
