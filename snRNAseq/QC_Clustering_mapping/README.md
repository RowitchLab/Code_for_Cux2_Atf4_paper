# QC, Clustering, and Cell Mapping Pipeline

This directory contains the quality control (QC), clustering, and cell type mapping scripts for the snRNA-seq analysis workflow. These scripts process raw 10X data through filtering, DESC clustering, and mapping to reference annotations.

## Overview

The pipeline follows this workflow:
1. **Data Loading** (`01_QC_desc_mapping.py`) - Load 10X h5 files and concatenate samples
2. **QC & Filtering** - Calculate QC metrics (mitochondrial %, ribosomal %, gene counts) and filter low-quality cells
3. **Doublet Detection** - Scrublet for doublet detection
4. **Normalization** - Log normalization and preprocessing
5. **DESC Clustering** - Deep Embedded Clustering for cell clustering
6. **Dimensionality Reduction** - UMAP generation at multiple resolutions
7. **Format Conversion** (`02_h5ad2rds.py`, `02_QC_R.py`) - Convert to Seurat format for R analysis

## Files

### Main Scripts

#### `01_QC_desc_mapping.py`
Complete QC and clustering pipeline for snRNA-seq data.

**Key Features:**
- Loads multiple 10X filtered h5 files and concatenates them
- Assigns sample metadata (Condition, Channel)
- Calculates QC metrics:
  - Mitochondrial gene percentage (`pct_counts_mt`)
  - Ribosomal gene percentage (`pct_counts_ribo`)
  - Hemoglobin gene percentage (`pct_counts_hb`)
- **Filtering thresholds:**
  - `min_genes=200` - minimum genes per cell
  - `pct_counts_mt < 10` - mitochondrial content filter
  - `min_cells=3` - minimum cells expressing a gene
- Doublet detection using Scrublet
- DESC clustering with multiple resolutions (0.6, 1.0, 1.2)
- UMAP visualization at each resolution

**Current Version Parameters (v03):**
```python
n_top_genes = 1024           # Highly variable genes
ae_dims = [1024, 256, 64]    # Auto-encoder dimensions
louvain_resolution = [0.6, 1.0, 1.2]  # Clustering resolutions
tneighbor = 20               # t-SNE/UMAP neighbor parameter
```

**Output Files:**
- `../Atf4/data/processed/10_load.h5ad` - Concatenated raw data
- `../Atf4/data/processed/11_preprocess.h5ad` - After QC filtering
- `../Atf4/data/processed/12_desc_v03.h5ad` - Final clustered object
- `results/v03/` - DESC clustering plots and UMAP visualizations

#### `02_h5ad2rds.py`
Converts AnnData h5ad objects to Seurat RDS format for R analysis.

**Features:**
- Handles raw count preservation (from `.raw` or `layers['counts']`)
- Cleans gene names (removes underscores, handles duplicates)
- Transfers cell metadata (obs)
- Optionally transfers dimensional reductions (UMAP, t-SNE, PCA)
- Optional log normalization (LogNormalize, scale.factor=10000)

**Usage:**
```python
from convert_h5ad_to_rds import convert_h5ad_to_rds

convert_h5ad_to_rds(
    h5ad_path='data/processed/13_desc_v03.h5ad',
    rds_path='data/processed/13_desc_v03.rds',
    save_coords=True,    # Transfer UMAP/t-SNE coordinates
    normalize=True       # Apply log normalization
)
```

#### `02_QC_R.py`
Comprehensive QC analysis with R/Seurat integration via rpy2.

**Classes and Methods:**

**`ComprehensiveSingleCellQC`**
- `convert_h5ad_to_rds()` - Convert AnnData to Seurat RDS with QC metrics
- `calculate_qc_metrics()` - Load and verify QC metrics
- `generate_qc_plots_tutorial_style()` - Create violin/scatter/density plots
- `detect_outliers_and_filter()` - Cell and gene filtering with visualization
- `generate_comprehensive_qc_report()` - Full QC workflow (deprecated)
- `generate_comprehensive_qc_report_v2()` - Advanced PBMC3k tutorial-compliant workflow

**Key QC Metrics Calculated:**
- `nFeature_RNA` - genes per cell
- `nCount_RNA` - UMI counts per cell
- `percent.mt` - mitochondrial gene %
- `percent.ribo` - ribosomal gene %
- `percent.hb` - hemoglobin gene %
- `log10GenesPerUMI` - complexity score

**Filtering Parameters (adjustable):**
- `min_features=200` - minimum genes
- `max_features=4000` - maximum genes
- `max_mt_percent=5.0` - max mitochondrial %
- `min_umis=500` - minimum UMI
- `max_umis=50000` - maximum UMI
- `min_complexity=0.8` - minimum log10(genes/UMI)

**Output Plots:**
- `qc_violin_plots_raw.pdf` - Distribution of QC metrics
- `qc_scatter_plots_raw.pdf` - Feature-feature relationships
- `qc_density_plots_raw.pdf` - Density distributions
- `filtering_comparison.pdf` - Before/after filtering comparison

**Example Usage:**
```python
from comprehensive_qc import ComprehensiveSingleCellQC

qc = ComprehensiveSingleCellQC()
results = qc.generate_comprehensive_qc_report_v2(
    h5ad_path='data/13_desc_v03.h5ad',
    output_dir='qc_results/',
    min_features=200,
    max_features=4000,
    max_mt_percent=5.0,
    n_variable_features=2000,
    n_pcs=10,
    clustering_resolution=0.5
)
```

## Data Samples

The pipeline processes data from the Cux2-Atf4 transgenic mouse study:

| Sample ID | Genotype | Condition |
|-----------|----------|-----------|
| 72_97, 75_99, 74_98 | Cux2-WT / Atf4-WT | Control WT |
| 25_88, 19_84, 11_81 | Cux2-Cre / Atf4-WT | Control Cre |
| 50_94, 9_93, 51_95 | Cux2-Cre / Atf4-flox | Cre/WT |
| 6_4, 10_92, 8_90 | Cux2-Cre / Atf4-flox | Cre/Cre |

## Dependencies

**Python Packages:**
```
scanpy
anndata
numpy
pandas
scipy
rpy2
```

**R Packages (via rpy2):**
```
Seurat
Matrix
dplyr
ggplot2
patchwork
viridis
cowplot
```

## How to Run

### Step 1: Data Loading and Clustering
```bash
# From snRNAseq/QC_Clustering_mapping/ directory
python 01_QC_desc_mapping.py
```
This generates DESC clustering and outputs `12_desc_v03.h5ad`

### Step 2: Convert to Seurat Format
```bash
python 02_h5ad2rds.py
```
This converts the h5ad file to RDS format for R analysis, creating `13_desc_v03.rds`

### Step 3: Comprehensive QC Analysis (Optional)
```bash
python 02_QC_R.py
```
This generates detailed QC plots and filtering statistics

## Output Structure

```
snRNAseq/
├── QC_Clustering_mapping/
│   ├── results/
│   │   └── v03/
│   │       ├── umap_desc_*.pdf       # UMAP plots at each resolution
│   │       └── ...clustering plots
│   └── qc_results/                    # Optional QC output
│
└── Atf4/data/processed/
    ├── 10_load.h5ad                  # Concatenated raw data
    ├── 11_preprocess.h5ad            # After QC filtering
    ├── 12_desc_v03.h5ad              # Clustered data
    └── 13_desc_v03.rds               # Seurat format
```

## Notes

- **QC Metrics**: All scripts follow standard snRNA-seq QC practices (removal of low-quality cells/genes, doublet detection)
- **DESC Clustering**: Deep Embedded Clustering provides robust clustering without pre-specification of cluster number
- **Multiple Resolutions**: Results are saved at resolutions 0.6, 1.0, and 1.2 for flexibility
- **Batch Correction**: Scrublet doublet detection is run per sample (`batch_key="Channel"`)

## Contact & Citation

For questions about this pipeline, refer to the main paper documentation.

