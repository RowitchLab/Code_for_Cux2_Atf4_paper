# Cux2 Analysis - Figure Generation Scripts

This directory contains scripts for generating publication figures from Cux2 transgenic mouse snRNA-seq data. The analysis integrates data from the Lucas and Popko public datasets with module scoring for stress response pathways.

## Overview

The Cux2 analysis focuses on:
1. **Differential Gene Expression (DGE) Analysis** - Identify DEGs between disease and control samples
2. **Module Scoring** - Score cells for stress response pathways (DDR, UPR, ISR, IFN)
3. **Figure Generation** - Create split violin plots and boxplots for publication
4. **Statistical Analysis** - Wilcoxon tests and effect size calculations

## Files

### Main Analysis Scripts

#### `Fig1_I_boxplot.py`
Performs differential expression analysis and generates boxplot figures for Figure 1.

**Purpose:** Quantifies the number of significantly altered genes (DEGs) per cell cluster between disease and control conditions, with optional filtering for specific gene sets (e.g., DNA Damage Response genes).

**Key Parameters:**
```python
fname = '../../h5ad_data/Lucas_org.h5ad'    # Input data
tcluster = 'anno3'                           # Cell type annotation column
fout_prefix = 'DEGbox_Lucas_withDDR'        # Output file prefix
num_cells_sampled = 200                      # Cells per cluster per iteration
tpval_cutoff = 0.01                         # Wilcoxon p-value threshold
tlog2fc = 0.5                               # Log2 fold-change threshold
sratio = 1.5                                # Size ratio for sampling
ddr_file = '../utils/DDR_list_final.csv'   # Optional: restrict to DDR genes
n_repeats = 50                              # Bootstrap iterations
```

**Workflow:**
1. Load h5ad file with Lucas dataset
2. For each cell cluster:
   - Subsample `num_cells_sampled` cells `n_repeats` times (with random state)
   - Perform Wilcoxon rank-sum test between disease and control
   - Count significant DEGs (p < 0.01, |log2FC| > 0.5)
   - Calculate distribution of DEG counts across iterations
3. If `ddr_file` provided:
   - Filter DEGs to only those in the DDR gene list
   - Count DDR-specific DEGs per iteration
4. Generate boxplot showing DEG count distribution across clusters

**Input Data Structure:**
- Must contain columns: `anno3` (cell types), `diagnosis` (disease status)
- Reference group: 'Control' (first in `tgroups`)
- Test group: 'MS' (second in `tgroups`)

**Output:**
```
figures/DEGbox_Lucas_withDDR_MS_vs_Control_boxplot2.pdf
```

**Functions:**
- `get_degs_per_cluster()` - Performs Wilcoxon test and extracts DEGs
- `intersection()` - Filters gene lists

**Example Usage:**
```python
# Customize parameters and run
python Fig1_I_boxplot.py
```

#### `Fig1_plot_Lucas_splitV.R`
Generates split violin plots for specific genes across conditions (Lucas dataset).

**Purpose:** Creates publication-quality split violin plots showing gene expression distributions split by disease condition (Control vs. MS).

**Input:**
- Lucas Seurat object: `./data/Lucas.rds`
- Cell type subset: EN-L2-3-A and EN-L2-3-B (cortical excitatory neurons)
- Disease variable: `diagnosis` (Control, MS)

**Genes Plotted (Figure 1):**
```r
genes_to_plot <- c('HDAC9', 'USP11', 'ACTB', 'GPX4', 'HSPA1A', 'ATM', 
                   'APEX1', 'XRCC6', 'PARP1', 'RAD23B', 'TP53BP1')
```

**Features:**
- Split by disease condition with custom colors: `#7A7A7A` (gray), `#DC9F41` (gold)
- Points jittered for visibility (pt.size=0.02, alpha=0.5)
- Minimal theme with black axis lines
- Individual PDF output for each gene

**Output Files:**
```
figures/TSplit_ViolinPlot_HDAC9.pdf
figures/TSplit_ViolinPlot_USP11.pdf
... (one per gene)
```

**Usage:**
```r
# Set working directory and source
setwd("/path/to/Cux2")
source("./Fig1_plot_Lucas_splitV.R")
```

#### `Fig3_popko_svlnplot.R`
Generates split violin plots for module scores across developmental stages (Popko dataset).

**Purpose:** Visualizes module scores for stress response pathways (DDR, UPR, ISR, IFN) across developmental ages in oligodendrocyte lineage cells.

**Input:**
- Popko Seurat object: `./data/Popko_EN-L2-3.rds`
- Cell subset: Ages W5/7, W17, W27/29, W41/44 (developmental timepoints)

**Modules Calculated:**
```r
# Module scores computed using AddModuleScore() with 500 control genes
seurat_subset <- add_module_score(seurat_subset, '../utils/DDR_list_final.csv', 'DDR')
seurat_subset <- add_module_score(seurat_subset, '../utils/IFN_sorted.csv', 'IFN')
seurat_subset <- add_module_score(seurat_subset, '../utils/UPR_list.csv', 'UPR')
seurat_subset <- add_module_score(seurat_subset, '../utils/ISR_list.csv', 'ISR')
```

**Gene Panels:**
1. **Module Scores** (initial analysis):
   - DDR1, UPR1, ISR1, IFN1

2. **Individual Genes** (main Figure 3):
   ```r
   genes_to_plot <- c('Atf4', 'Cux2', 'Hdac9', 'Actb', 'Gpx4', 
                      '2900097C17Rik', 'Apex1', 'Rad23b', 'Trp53bp1', 'Xrcc6')
   ```

**Features:**
- **Split violin plots** separating conditions across age groups
- **Statistical testing** via `summarize_group_stats()` (Wilcoxon rank-sum test)
- **Custom colors** via `my_palette()` function
- **Output formats:**
  - PDF split violin plots: `figures/TSplit_ViolinPlot_[GENE]_Popko_ENL2-3.pdf`
  - CSV statistics: `csv/[GENE]_group_stats_Popko_ENL2-3.csv`

**Statistical Output:**
Each gene generates a CSV with:
- Mean expression per condition/timepoint
- Cell counts per group
- Wilcoxon p-values (formatted as scientific notation)

**Usage:**
```r
# From Cux2 directory
setwd("/path/to/Cux2")
source("./Fig3_popko_svlnplot.R")
```

## Utility Scripts

### `../utils/splitviolin_plot_base.R`
Base functions for split violin plot generation and statistics.

**Key Functions:**
- `add_module_score()` - Adds module scores to Seurat object using gene lists
- `plot_split_violin()` - Generates split violin plots with Wilcoxon tests
- `summarize_group_stats()` - Calculates and saves group statistics (mean expression, cell counts, p-values)
- `my_palette()` - Custom color palette for plots

**Features:**
- Uses `ggplot2::geom_split_violin()` for split violin visualization
- Incorporates `ggpubr::stat_compare_means()` for p-value annotation
- Automatic statistical testing (Wilcoxon rank-sum test)

### `../utils/DDR_list_final.csv`
Comprehensive DNA Damage Response gene list (~500 genes).

**Contents:**
Genes involved in DNA repair pathways including:
- Homologous recombination (RAD51, BRCA1, BRCA2, etc.)
- Non-homologous end joining (XRCC4, XRCC5, LIG4, etc.)
- Base excision repair (APEX1, OGG1, UNG, etc.)
- Nucleotide excision repair (ERCC1-6, GTF2H genes, etc.)
- Mismatch repair (MLH1, MLH3, MSH2-6, etc.)
- And many more DNA damage response genes

**Usage in Scripts:**
```python
ddr_file = '../utils/DDR_list_final.csv'
ddr_list = pd.read_csv(ddr_file, header=None)[0].tolist()
```

### `../utils/box_plot_base.R`
Base function for generating boxplots with sorted categories.

**Function:** `generate_boxplot()`
- Reshapes data from wide to long format
- Sorts clusters by mean value
- Creates notched boxplot with outliers shown
- Saves as PDF with customizable dimensions

## Data Files

### Lucas Dataset
- **File:** `./data/Lucas.rds`
- **Cell types:** EN-L2-3-A, EN-L2-3-B (cortical excitatory neurons)
- **Conditions:** Control, MS (Multiple Sclerosis)
- **Metadata columns:** `cell_type`, `diagnosis`

### Popko Dataset
- **File:** `./data/Popko_EN-L2-3.rds`
- **Cell types:** Oligodendrocyte lineage cells (EN-L2-3 region)
- **Ages:** W5/7, W17, W27/29, W41/44 (weeks post-natal)
- **Metadata columns:** `age`, `Condition`, and module scores

## Workflow Summary

```
Figure Generation Pipeline:
├── Fig1_I_boxplot.py
│   ├── Load Lucas h5ad data
│   ├── Calculate DEG counts per cluster (50 iterations)
│   ├── Optional: Filter to DDR genes
│   └── Output: Boxplot PDF
│
├── Fig1_plot_Lucas_splitV.R
│   ├── Load Lucas Seurat object
│   ├── Subset EN-L2-3 cells
│   ├── Generate split violin plots (11 genes)
│   └── Output: 11 PDF files
│
└── Fig3_popko_svlnplot.R
    ├── Load Popko Seurat object
    ├── Subset developmental timepoints
    ├── Calculate module scores
    ├── Generate split violin plots (10 genes)
    ├── Calculate statistics
    └── Output: 10 PDFs + 10 CSVs
```

## Dependencies

**Python Packages:**
- scanpy
- pandas
- numpy
- rpy2 (for R integration)

**R Packages:**
- Seurat
- ggplot2
- introdataviz (for `geom_split_violin()`)
- forcats
- ggpubr

## How to Run

### Figure 1 Analysis (Python & R)
```bash
# From Cux2 directory
# Step 1: Generate DEG boxplot (Python)
python Fig1_I_boxplot.py

# Step 2: Generate split violin plots (R)
Rscript Fig1_plot_Lucas_splitV.R
```

### Figure 3 Analysis (R)
```bash
# From Cux2 directory
Rscript Fig3_popko_svlnplot.R
```

## Output Organization

```
Cux2/
├── figures/
│   ├── TSplit_ViolinPlot_HDAC9.pdf
│   ├── TSplit_ViolinPlot_USP11.pdf
│   ├── ... (Lucas split violin plots)
│   ├── TSplit_ViolinPlot_Atf4_Popko_ENL2-3.pdf
│   ├── ... (Popko split violin plots)
│   └── DEGbox_Lucas_withDDR_MS_vs_Control_boxplot2.pdf
│
└── csv/
    ├── Atf4_group_stats_Popko_ENL2-3.csv
    ├── Cux2_group_stats_Popko_ENL2-3.csv
    └── ... (statistics for each gene)
```

## Notes

- **Sampling Strategy:** Fig1_I_boxplot uses bootstrap resampling (n=50) to estimate variability
- **Statistical Tests:** All comparisons use non-parametric Wilcoxon rank-sum test (appropriate for non-normal distributions)
- **Module Scoring:** Uses Seurat's `AddModuleScore()` with 500 control genes for robust scoring
- **Color Scheme:** Consistent across figures (gray/gold for controls/conditions)

## Citation & Related Analysis

These figures are part of the Cux2/ATF4 transgenic mouse study.
Related analyses in:
- QC_Clustering_mapping/ - Data preprocessing and clustering
- Atf4/ - Public dataset validation

