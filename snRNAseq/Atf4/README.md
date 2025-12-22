# ATF4 Analysis - Public Dataset Integration

This directory contains analysis scripts for exploring ATF4 expression in the Linnarsson public datasets, with integration of developmental and adolescent brain transcriptomics data.

## Overview

This analysis focuses on:
1. Loading and converting Linnarsson public loom files to AnnData format
2. Subsetting cell types and developmental stages
3. Converting to Seurat format for R analysis
4. Generating publication-ready violin plots showing ATF4 expression patterns across cell types and developmental ages

## Files

### Data Subdirectory

#### `public_dataset/01_load_convert.py`
Loads and converts Linnarsson public dataset from loom to h5ad format.

**Datasets Processed:**
1. **Linnarsson Adolescent Data** (`linnarsson_adolescent.loom`)
   - Multiple ages: p21, p22, p23, p25, p26, p27, p28, p29, p60
   - Cell types: Telencephalon inhibitory interneurons, projecting excitatory neurons
   
2. **Linnarsson Developmental Data** (`linnarsson_dev_all.loom`)
   - Embryonic ages: e17.0, e17.5, e18.0
   - Cell type subclasses:
     - Cerebellar glutamatergic/GABAergic
     - Forebrain glutamatergic/GABAergic
     - Midbrain glutamatergic/GABAergic
     - Cortical/hippocampal glutamatergic

**Key Processing Steps:**
- Extracts UMAP coordinates from metadata (fields: `_X`, `_Y`)
- Extracts t-SNE coordinates (`_tSNE1`, `_tSNE2`)
- Extracts PCA coordinates (`_PC1`, `_PC2`)
- Preserves counts matrix in `layers["counts"]`
- Normalizes and log-transforms data
- Outputs standardized h5ad format

**Output Files:**
- `../h5ad_data/linnarsson_adolescent.h5ad` (~1.2 GB)
- `../h5ad_data/linnarsson_dev_all.h5ad` (~1.5 GB)

**Usage:**
```bash
cd public_dataset/
python 01_load_convert.py
```

#### `public_dataset/load_data.py`
Loads processed h5ad files and converts to Seurat RDS format with additional processing.

**Key Functions:**
- `convert_h5ad_to_rds()` - Converts AnnData to Seurat RDS
  - Preserves raw count matrices
  - Cleans gene names (removes underscores, handles duplicates)
  - Transfers UMAP/t-SNE embeddings as Seurat reductions
  - Optional log normalization

**QC Metrics Added:**
- Mitochondrial gene % (`percent.mt`)
- Ribosomal gene % (`percent.ribo`)
- Hemoglobin gene % (`percent.hb`)
- Platelet markers % (`percent.platelet`)
- Gene complexity score (`log10GenesPerUMI`)

**Output Files:**
- `linnarsson_adolescent_filtered.rds` - Processed adolescent data
- `linnarsson_dev_all_filtered.rds` - Processed developmental data

#### `public_dataset/load_data.R`
R script for loading Linnarsson loom files directly using Seurat/SeuratDisk.

**Requirements:**
```r
library(Seurat)
library(SeuratDisk)
library(loomR)
```

#### `public_dataset/plot_violin.R`
Generates publication-quality violin plots of ATF4 expression patterns.

**Features:**
- Uses `plot1cell::complex_vlnplot_single()` for advanced plotting
- Creates split violin plots showing ATF4 across:
  - Developmental ages
  - Cell type subclasses
- High-quality PDF output with customizable themes

**Key Plots Generated:**
- `violin_Age_Atf4.pdf` - ATF4 by age in adolescent data
- `violin_Age_Atf4_dev.pdf` - ATF4 by cell type subclass across developmental stages
- `violin_Age_Atf4_dev2.pdf` - Alternative layout with rotated axis labels

**Usage:**
```r
# Set working directory
setwd("./public_dataset")

# Load Seurat object
seurat_object <- readRDS("./linnarsson_adolescent_filtered.rds")

# Generate plot
pdf("violin_Age_Atf4.pdf", width=5, height=4)
complex_vlnplot_single(seurat_object, feature = "Atf4", groups = "Age")
dev.off()
```

### Analysis Workflow

```
Data Processing Pipeline:
├── 01_load_convert.py
│   ├── Load linnarsson_*.loom
│   ├── Extract coordinates (UMAP, t-SNE, PCA)
│   ├── Normalize and log-transform
│   └── Output: h5ad files
│
├── load_data.py / load_data.R
│   ├── Load h5ad files
│   ├── Calculate QC metrics
│   ├── Convert to Seurat RDS
│   └── Output: .rds files
│
└── plot_violin.R
    ├── Load Seurat objects
    ├── Filter by age/subclass
    ├── Generate violin plots
    └── Output: PDF figures
```

## Datasets

### Linnarsson Adolescent Data
**Source:** Linnarsson lab, adolescent mouse brain
- **Ages:** Postnatal days 21-60
- **Cell types:** Telencephalon interneurons and excitatory neurons
- **Total cells:** ~100,000+
- **Key metadata:** Age, TaxonomyRank4 (cell type)

### Linnarsson Developmental Data
**Source:** Linnarsson lab, embryonic and early postnatal mouse brain
- **Ages:** Embryonic days E17.0-E18.0
- **Cell types:** Multiple subclasses across brain regions (cerebellum, forebrain, midbrain)
- **Total cells:** ~200,000+
- **Key metadata:** Age, Subclass, Class

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

**R Packages:**
```
Seurat
SeuratDisk
loomR
plot1cell
ggplot2
```

## How to Run

### Full Pipeline (Python)
```bash
# Change to public_dataset directory
cd public_dataset/

# Step 1: Load and convert loom files to h5ad
python 01_load_convert.py

# Step 2: Convert to Seurat RDS format
python load_data.py
```

### Plotting (R)
```r
# Set working directory
setwd("./public_dataset")
source("plot_violin.R")

# Plots will be generated as PDF files
```

## Output Files

```
Atf4/
├── data/
│   └── h5ad_data/
│       ├── linnarsson_adolescent.h5ad        # ~1.2 GB
│       ├── linnarsson_adolescent.h5ad        # ~1.5 GB
│       ├── linnarsson_adolescent_filtered.rds # Converted
│       └── linnarsson_dev_all_filtered.rds   # Converted
│
└── public_dataset/
    ├── violin_Age_Atf4.pdf                  # Main figure
    ├── violin_Age_Atf4_dev.pdf               # Developmental figure
    ├── violin_Age_Atf4_dev2.pdf              # Alternative layout
    └── linnarsson_*_filtered.rds
```

## Filtering & Subsetting

### Adolescent Data Subsetting
```python
adata = adata[adata.obs.Age.isin(["p21", "p22", "p23", "p25", "p26", "p27", "p28", "p29", "p60"])]
adata = adata[adata.obs['TaxonomyRank4'].isin(['Telencephalon inhibitory interneurons',
                                                'Telencephalon projecting excitatory neurons'])]
```

### Developmental Data Subsetting
```python
adata = adata[adata.obs.Age.isin(["e17.0", "e18.0", "e17.5"])]
adata = adata[adata.obs['Subclass'].isin(["Cerebelllum glutamatergic",
                                          "Cerebelllum GABAergic",
                                          "Forebrain glutamatergic",
                                          "Forebrain GABAergic",
                                          "Midbrain glutamatergic",
                                          "Midbrain GABAergic",
                                          "Cortical or hippocampal glutamatergic"])]
```

## QC & Quality Notes

- All data is subset to specific cell types of interest before conversion
- Raw count matrices are preserved for reproducibility
- Normalized data is log-transformed for visualization
- Dimensional reductions (UMAP, t-SNE, PCA) are transferred to Seurat format

## Related Analysis

This public dataset analysis provides:
- **Developmental context** for ATF4 expression patterns
- **Cell-type specificity** of ATF4 induction
- **Age-dependent trends** in ATF4 expression
- **Comparison baseline** for transgenic mouse data

## References

**Linnarsson Lab Datasets:**
- Adolescent: https://www.linnarsson.org/blog/single-cell-rna-seq-databases
- Developmental: Available through the same resource

## Contact

For questions about data processing or figure generation, refer to the main analysis documentation.

