# snRNA-seq Analysis Pipeline - Cux2/ATF4 Study

Comprehensive single-nucleus RNA-sequencing analysis workflow for the Cux2/ATF4 transgenic mouse study. This pipeline covers data preprocessing, quality control, clustering, and figure generation for publication.

## Project Overview

This analysis investigates the role of ATF4 (Activating Transcription Factor 4) in Cux2+ cortical excitatory neurons during development and in disease conditions (Multiple Sclerosis model). The study integrates:
- **Transgenic mouse data** - Cux2-Cre driven ATF4 conditional knockout
- **Public datasets** - Linnarsson developmental and adolescent brain data
- **Published datasets** - Lucas (MS) and Popko (developmental oligodendrocytes) data

## Repository Structure

```
snRNAseq/
│
├── QC_Clustering_mapping/          # Data preprocessing & clustering
│   ├── README.md                   # Detailed documentation
│   ├── 01_QC_desc_mapping.py      # Load, QC filter, DESC clustering
│   ├── 02_h5ad2rds.py              # AnnData → Seurat conversion
│   ├── 02_QC_R.py                 # Comprehensive R-based QC analysis
│   └── results/                    # Clustering outputs (UMAP plots)
│
├── Atf4/                           # ATF4-specific analyses
│   ├── README.md                   # Documentation
│   ├── data/                       # Raw and processed data
│   │   ├── processed/
│   │   │   ├── 10_load.h5ad       # Concatenated raw data
│   │   │   ├── 11_preprocess.h5ad # After QC filtering
│   │   │   ├── 12_desc_v03.h5ad   # Final clustered data
│   │   │   └── 13_desc_v03.rds    # Seurat format
│   │   └── h5ad_data/              # Public datasets
│   │       ├── linnarsson_adolescent.h5ad
│   │       └── linnarsson_dev_all.h5ad
│   │
│   └── public_dataset/             # Public data integration
│       ├── 01_load_convert.py      # Loom → h5ad conversion
│       ├── load_data.py            # h5ad → Seurat RDS
│       ├── load_data.R             # R-based loading
│       └── plot_violin.R           # ATF4 expression plots
│
├── Cux2/                           # Cux2 transgenic analysis & figures
│   ├── README.md                   # Documentation
│   ├── data/                       # Seurat objects
│   │   ├── Lucas.rds              # Lucas MS dataset
│   │   └── Popko_EN-L2-3.rds      # Popko developmental dataset
│   │
│   ├── Fig1_I_boxplot.py           # DEG boxplot (Figure 1I)
│   ├── Fig1_plot_Lucas_splitV.R   # Split violin plots (Figure 1)
│   ├── Fig3_popko_svlnplot.R      # Module score plots (Figure 3)
│   │
│   ├── figures/                    # Output plots (PDF)
│   └── csv/                        # Statistical outputs
│
├── utils/                          # Shared utilities
│   ├── DDR_list_final.csv         # DNA Damage Response gene list
│   ├── get_go_lists.py             # Gene list utilities
│   ├── splitviolin_plot_base.R    # Base plot functions
│   ├── box_plot_base.R             # Boxplot utilities
│   └── *.csv                       # Other gene lists (IFN, UPR, ISR)
│
└── ReadMe.md                       # This file
```

## Analysis Pipeline

### Phase 1: Data Preprocessing and Clustering
**Location:** `QC_Clustering_mapping/`

```
Raw 10X h5 files (6 samples)
         ↓
01_QC_desc_mapping.py
  ├─ Load & concatenate 10X files
  ├─ Calculate QC metrics (mitochondrial %, genes/cell, UMI/cell)
  ├─ Filter cells: min_genes=200, mt%<10, min_cells=3 genes
  ├─ Doublet detection (Scrublet)
  ├─ Normalization & log-transform
  └─ DESC clustering (resolutions: 0.6, 1.0, 1.2)
         ↓
11_preprocess.h5ad (filtered, normalized)
         ↓
02_h5ad2rds.py + 02_QC_R.py
  ├─ Convert to Seurat RDS format
  ├─ Generate comprehensive QC plots
  └─ Transfer dimensional reductions (UMAP)
         ↓
13_desc_v03.rds (Seurat object, ready for analysis)
```

**Outputs:**
- Clustering results at 3 resolutions (0.6, 1.0, 1.2)
- UMAP visualizations
- QC metrics and filtering statistics
- Seurat-compatible RDS file for downstream analysis

### Phase 2: Public Data Integration
**Location:** `Atf4/public_dataset/`

```
Linnarsson Loom files (adolescent + developmental)
         ↓
01_load_convert.py
  ├─ Extract metadata & embeddings
  ├─ Subset to cell types of interest
  ├─ Normalize & log-transform
  └─ Output standardized h5ad
         ↓
linnarsson_*.h5ad
         ↓
load_data.py / load_data.R
  ├─ Calculate QC metrics (Seurat)
  ├─ Transfer embeddings
  └─ Convert to RDS
         ↓
linnarsson_*_filtered.rds (Ready for plotting)
         ↓
plot_violin.R
  └─ Generate ATF4 expression plots
```

**Outputs:**
- Standardized AnnData and Seurat objects from public data
- Publication-quality violin plots of ATF4 expression
- Developmental and cell-type-specific expression patterns

### Phase 3: Figure Generation
**Location:** `Cux2/`

#### Figure 1: Lucas MS Dataset
```
Lucas.rds (EN-L2-3 cells, Control vs MS)
         ↓
Split Violin Plots (Fig1_plot_Lucas_splitV.R)
  └─ 11 genes: HDAC9, USP11, ACTB, GPX4, HSPA1A, ATM, 
               APEX1, XRCC6, PARP1, RAD23B, TP53BP1

DEG Boxplot (Fig1_I_boxplot.py)
  ├─ Bootstrap DEG analysis (50 iterations, n=200 cells/cluster)
  ├─ Wilcoxon test (p<0.01, |log2FC|>0.5)
  └─ Optional: Filter to DDR genes
```

**Outputs:**
- 11 individual split violin plot PDFs
- Boxplot showing DEG count variability across clusters

#### Figure 3: Popko Developmental Data
```
Popko_EN-L2-3.rds (4 developmental timepoints)
         ↓
Module Score Calculation
  ├─ DDR: DNA Damage Response
  ├─ UPR: Unfolded Protein Response
  ├─ ISR: Integrated Stress Response
  └─ IFN: Interferon response

Split Violin Plots (Fig3_popko_svlnplot.R)
  ├─ Module scores (4 features)
  ├─ Individual genes (10 features)
  │   └─ Atf4, Cux2, Hdac9, Actb, Gpx4, 2900097C17Rik,
  │      Apex1, Rad23b, Trp53bp1, Xrcc6
  └─ Statistical testing (Wilcoxon, with p-values)
```

**Outputs:**
- Split violin plots (PDF) for module scores and genes
- CSV tables with statistics (mean expression, p-values, cell counts)

## Data Summary

| Dataset | Cell Types | Conditions | N Samples | N Cells | Source |
|---------|-----------|-----------|-----------|---------|--------|
| Transgenic Cux2/ATF4 | Multiple | 4 genotypes | 6 | ~100K | Generated |
| Linnarsson Adolescent | Telencephalon neurons | Development | - | ~100K+ | Public |
| Linnarsson Developmental | Multi-region | E17-E18 | - | ~200K+ | Public |
| Lucas | Cortical neurons (EN-L2-3) | Control / MS | - | ~50K | Public |
| Popko | Oligodendrocytes | Development | - | ~50K | Public |

## Key Parameters

### QC & Filtering
```
min_genes = 200               # Minimum genes/cell
max_genes = 4000 (optional)   # Maximum genes/cell
min_cells = 3                 # Minimum cells/gene
pct_mt < 10%                  # Mitochondrial gene filter
```

### Clustering (DESC)
```
n_top_genes = 1024            # Highly variable genes
ae_dims = [1024, 256, 64]     # Auto-encoder dimensions
resolutions = [0.6, 1.0, 1.2] # Louvain resolutions
n_neighbors = 20              # t-SNE/UMAP neighbor parameter
```

### Figure Generation
```
Sampling: n_repeats = 50, n_cells = 200/cluster
Statistical test: Wilcoxon rank-sum (non-parametric)
Fold-change threshold: |log2FC| > 0.5
P-value threshold: p < 0.01
Module scoring: AddModuleScore with 500 control genes
```

## Required Software & Packages

### Python Environment
```
Python 3.8+
scanpy >= 1.8
anndata >= 0.8
pandas >= 1.1
numpy >= 1.19
scipy >= 1.5
scikit-learn >= 0.24
rpy2 >= 3.4          # For R integration
```

### R Environment
```
R >= 4.0
Seurat >= 4.0
ggplot2 >= 3.3
dplyr >= 1.0
patchwork >= 1.1
introdataviz         # For split violin plots
ggpubr >= 0.4        # For statistical annotations
loomR                 # For reading loom files
SeuratDisk           # For H5Seurat conversion
```

## Installation & Setup

### Python
```bash
# Create conda environment
conda create -n cux2_atf4 python=3.9 -y
conda activate cux2_atf4

# Install packages
pip install scanpy anndata pandas scipy rpy2
```

### R
```r
# Install core packages
install.packages(c("ggplot2", "dplyr", "patchwork", "ggpubr"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Seurat")

# Install introdataviz (for geom_split_violin)
devtools::install_github("tejasjinturkar/introdataviz")
```

## How to Run - Complete Workflow

### 1. Preprocessing & Clustering
```bash
cd QC_Clustering_mapping/
python 01_QC_desc_mapping.py        # ~1-2 hours
python 02_h5ad2rds.py               # ~15 minutes
```

### 2. Public Data Integration
```bash
cd ../Atf4/public_dataset/
python 01_load_convert.py            # ~30 minutes (first time, downloads data)
python load_data.py                  # ~15 minutes
Rscript plot_violin.R                # ~5 minutes
```

### 3. Figure Generation
```bash
cd ../../Cux2/
python Fig1_I_boxplot.py             # ~30 minutes
Rscript Fig1_plot_Lucas_splitV.R    # ~5 minutes
Rscript Fig3_popko_svlnplot.R       # ~5 minutes
```

## Output Files & Locations

### Data Objects
- `Atf4/data/processed/12_desc_v03.h5ad` - Final clustered AnnData (1-2 GB)
- `Atf4/data/processed/13_desc_v03.rds` - Seurat format (~800 MB)
- `Atf4/data/h5ad_data/*.h5ad` - Public datasets (1-2 GB each)

### Figures
- `QC_Clustering_mapping/results/v03/` - UMAP clustering plots
- `Atf4/public_dataset/violin_*.pdf` - ATF4 expression plots
- `Cux2/figures/` - Publication-ready figures (split violins, boxplots)

### Statistics
- `Cux2/csv/` - Gene expression statistics and p-values
- `QC_Clustering_mapping/qc_results/` - QC metrics and filtering statistics

## Troubleshooting

### Common Issues

**Issue:** "R_HOME not found" error
```bash
# Set R environment variables
export R_HOME=/path/to/R
export R_LD_LIBRARY_PATH=/path/to/R/lib
```

**Issue:** rpy2 conversion errors
- Ensure R packages are installed and compatible with rpy2 version
- Reinstall rpy2: `pip install --upgrade rpy2`

**Issue:** Loom file too large
- Subsetting in 01_load_convert.py can reduce memory usage
- Process in chunks if needed

## Citation

Please cite the following when using this pipeline:
- Original analysis: [Citation here]
- Seurat: Stuart et al., Cell (2019)
- scanpy: Wolf et al., Genome Biology (2018)

## Contact

For questions or issues with the analysis pipeline, please contact [lab/author information].

## Notes

- **Reproducibility:** All random seeds are set for reproducibility. Use `random_state` parameters in Python and `set.seed()` in R
- **Computational Requirements:** 
  - Total runtime: ~2-3 hours (depending on machine)
  - Memory: 16+ GB RAM recommended
  - Disk space: ~10 GB for raw + processed data
- **Version Control:** Analysis conducted with desc_pytorch package for clustering
- **Data Management:** Large h5ad/rds files should be stored separately and referenced via paths

