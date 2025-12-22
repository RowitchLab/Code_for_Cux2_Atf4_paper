import os
os.environ["R_HOME"]="/Library/Frameworks/R.framework/Resources"
os.environ["R_LD_LIBRARY_PATH"]="/Library/Frameworks/R.framework/Resources/lib"
os.environ["RSCRIPT_PATH"]="/usr/local/bin/Rscript"

import anndata
import scipy.sparse
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import StrVector, FloatMatrix, ListVector, FloatVector
from typing import Dict, List, Tuple, Optional, Union
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

class ComprehensiveSingleCellQC:
    """
    Complete Single Cell QC workflow using Seurat through rpy2
    Based on https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
    """
    
    def __init__(self):
        """Initialize R packages for QC analysis"""
        try:
            # Import R packages
            self.base = importr('base')
            self.seurat = importr('Seurat')
            self.Matrix = importr('Matrix')
            self.dplyr = importr('dplyr')
            self.ggplot2 = importr('ggplot2')
            self.patchwork = importr('patchwork')
            self.viridis = importr('viridis')
            self.cowplot = importr('cowplot')
            print("✅ All required R packages loaded successfully")
            
        except Exception as e:
            print(f"❌ Error loading R packages: {e}")
            print("Please install required R packages:")
            print("install.packages(c('Seurat', 'dplyr', 'ggplot2', 'patchwork', 'viridis', 'cowplot'))")
    
    def convert_h5ad_to_rds(self, h5ad_path: str, rds_path: str, save_coords=False, normalize=False):
        """
        Convert AnnData h5ad file to Seurat RDS format
        Based on your provided conversion function
        """
        print(f"🔄 Converting {h5ad_path} to Seurat format...")
        
        # Read h5ad file
        adata = anndata.read_h5ad(h5ad_path)
        print(f'Original shape: {adata.shape}')

        # Get raw counts matrix
        if adata.raw is not None:
            counts = adata.raw.X.copy()
            var_names = adata.raw.var_names
            print('Using adata.raw.X as raw counts')
        elif 'counts' in adata.layers:
            counts = adata.layers['counts'].copy()
            var_names = adata.var_names
            print('Using adata.layers["counts"] as raw counts')
        else:
            counts = adata.X.copy()
            var_names = adata.var_names
            print('WARNING: no raw counts found -> using adata.X')

        # Ensure sparse matrix
        if not scipy.sparse.issparse(counts):
            counts = scipy.sparse.csr_matrix(counts)

        # Remove duplicate gene names
        var_names = pd.Index(var_names).astype(str)
        if var_names.duplicated().any():
            var_names = var_names.make_unique()
        counts = counts[:, ~var_names.duplicated()].copy()
        adata = adata[:, var_names].copy()
        print(f'After duplicate removal: {adata.shape}')

        # Create R sparse matrix (genes x cells for Seurat)
        counts_coo = counts.T.tocoo()
        genes = var_names.tolist()
        cells = adata.obs_names.astype(str).tolist()

        r_counts = self.Matrix.sparseMatrix(
            i=ro.IntVector(counts_coo.row + 1),
            j=ro.IntVector(counts_coo.col + 1),
            x=ro.FloatVector(counts_coo.data),
            dims=ro.IntVector(counts_coo.shape),
            dimnames=ro.r.list(genes=StrVector(genes), cells=StrVector(cells))
        )

        # Clean gene names (remove underscores, ensure uniqueness)
        ro.r.assign('counts', r_counts)
        ro.r('''
        rn <- rownames(counts)
        rn[rn == ""] <- NA
        rn <- gsub("_", "-", rn)
        rn[is.na(rn)] <- paste0("GENE-", seq_len(sum(is.na(rn))))
        rownames(counts) <- make.unique(rn)
        ''')
        r_counts_clean = ro.r('counts')

        # Convert metadata
        with localconverter(default_converter + pandas2ri.converter):
            meta_r = pandas2ri.py2rpy(adata.obs)

        # Create Seurat object
        ro.r.assign('counts_clean', r_counts_clean)
        ro.r.assign('meta_data', meta_r)
        ro.r('seurat_obj <- CreateSeuratObject(counts = counts_clean, meta.data = meta_data)')
        seurat_obj = ro.r('seurat_obj')

        # Add dimensional reductions if requested
        if save_coords:
            ro.r.assign('seurat_obj', seurat_obj)
            ro.r('if (!"reductions" %in% slotNames(seurat_obj)) seurat_obj@reductions <- list()')
            
            embed_map = {'X_umap': 'UMAP', 'X_tsne': 'tSNE', 'X_pca': 'PCA'}
            for key, seurat_name in embed_map.items():
                if key in adata.obsm:
                    coord = adata.obsm[key][:, :2] if adata.obsm[key].shape[1] >= 2 else adata.obsm[key]
                    coord_r = ro.r.matrix(FloatVector(coord.T.ravel()),
                                        nrow=coord.shape[0],
                                        ncol=coord.shape[1],
                                        dimnames=ro.r.list(StrVector(cells),
                                                          StrVector([f'Dim{i+1}' for i in range(coord.shape[1])])))
                    
                    ro.r.assign('seurat_obj', seurat_obj)
                    ro.r.assign('coord_r', coord_r)
                    ro.r.assign('name', seurat_name)
                    ro.r('seurat_obj@reductions[[name]] <- CreateDimReducObject(embeddings = coord_r, key = paste0(name, "_"), assay = DefaultAssay(seurat_obj))')
                    seurat_obj = ro.r('seurat_obj')  # Get updated object
                    print(f'  ✅ Added {seurat_name} reduction')

        # Normalize data if requested
        if normalize:
            ro.r.assign('seurat_obj', seurat_obj)
            ro.r('seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)')
            print('✅ Data normalized (LogNormalize)')
            seurat_obj = ro.r('seurat_obj')

        # Calculate QC metrics immediately and save them to the object
        print("🔬 Calculating and adding QC metrics to Seurat object...")
        ro.r.assign('seurat_obj', seurat_obj)
        
        qc_calc_code = '''
        library(Seurat)
        
        # Calculate mitochondrial gene percentage
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
        
        # Calculate ribosomal gene percentage  
        seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
        
        # Calculate hemoglobin gene percentage
        seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(p)]")
        
        # Calculate platelet gene percentage
        seurat_obj[["percent.platelet"]] <- PercentageFeatureSet(seurat_obj, pattern = "PPBP|PF4|SDPR|GP1BA|TUBB1|CLU|GNG11|NRGN|RGS18|TPM4")
        
        # Calculate complexity score
        seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
        
        # Calculate mitochondrial ratio
        seurat_obj$mitoRatio <- seurat_obj@meta.data$percent.mt / 100
        
        # Print summary of calculated metrics
        cat("QC metrics calculated:\\n")
        cat("  Percent MT: mean =", round(mean(seurat_obj$percent.mt), 2), "%\\n")
        cat("  Percent Ribo: mean =", round(mean(seurat_obj$percent.ribo), 2), "%\\n") 
        cat("  Percent HB: mean =", round(mean(seurat_obj$percent.hb), 2), "%\\n")
        cat("  Complexity: mean =", round(mean(seurat_obj$log10GenesPerUMI, na.rm = TRUE), 3), "\\n")
        '''
        
        ro.r(qc_calc_code)
        seurat_obj = ro.r('seurat_obj')
        print('✅ QC metrics calculated and added to object')

        # Save RDS file with QC metrics included
        ro.r.assign('seurat_obj_final', seurat_obj)
        ro.r(f'saveRDS(seurat_obj_final, file = "{rds_path}")')
        print(f'✅ Seurat object with QC metrics saved to: {rds_path}')
        
        return seurat_obj, adata
    
    def calculate_qc_metrics(self, seurat_obj_path: str):
        """
        Calculate comprehensive QC metrics following the tutorial approach
        This function now just returns the metadata since QC metrics are calculated during conversion
        """
        print("📊 Loading QC metrics from Seurat object...")
        
        # Load Seurat object
        ro.r(f'seurat_obj <- readRDS("{seurat_obj_path}")')
        
        # Verify QC metrics are present
        check_code = '''
        available_cols <- colnames(seurat_obj@meta.data)
        cat("Available metadata columns:\\n")
        cat(paste(available_cols, collapse = ", "), "\\n\\n")
        
        # Check specific QC metrics
        qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "log10GenesPerUMI")
        for (metric in qc_metrics) {
            if (metric %in% available_cols) {
                cat("✅", metric, "- available\\n")
            } else {
                cat("❌", metric, "- missing\\n")
            }
        }
        '''
        
        ro.r(check_code)
        print("✅ QC metrics verification complete")
        
        # Get metadata for analysis
        metadata_r = ro.r('seurat_obj@meta.data')
        with localconverter(default_converter + pandas2ri.converter):
            metadata = pandas2ri.rpy2py(metadata_r)
        
        return metadata
    
    def generate_qc_plots_tutorial_style(self, seurat_obj_path: str, output_dir: str):
        """
        Generate QC plots following the tutorial style
        """
        print("📊 Generating comprehensive QC plots...")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Load the Seurat object and ensure QC metrics are calculated
        ro.r(f'seurat_obj <- readRDS("{seurat_obj_path}")')
        
        # First, ensure all QC metrics are calculated
        qc_prep_code = '''
        library(Seurat)
        library(ggplot2)
        library(patchwork)
        
        # Make sure QC metrics are calculated
        if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
            seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
        }
        if (!"percent.ribo" %in% colnames(seurat_obj@meta.data)) {
            seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
        }
        if (!"log10GenesPerUMI" %in% colnames(seurat_obj@meta.data)) {
            seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
        }
        
        # Check what columns we have
        available_features <- colnames(seurat_obj@meta.data)
        cat("Available features:", paste(available_features, collapse = ", "), "\\n")
        '''
        
        ro.r(qc_prep_code)
        
        # Generate violin plots for basic metrics (always available)
        try:
            violin_plot_code = f'''
            # Violin plot for number of genes detected
            p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0.1) +
                  theme(legend.position = "none") +
                  labs(title = "Genes Detected per Cell", x = "", y = "Number of Genes")
            
            # Violin plot for total UMI counts
            p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1) +
                  theme(legend.position = "none") +
                  labs(title = "UMI Counts per Cell", x = "", y = "UMI Counts")
            
            # Check if percent.mt exists before plotting
            if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {{
                p3 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0.1) +
                      theme(legend.position = "none") +
                      labs(title = "Mitochondrial Gene %", x = "", y = "Percent MT")
            }} else {{
                p3 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No MT genes found") +
                      labs(title = "Mitochondrial Gene %") +
                      theme_minimal()
            }}
            
            # Check if percent.ribo exists before plotting
            if ("percent.ribo" %in% colnames(seurat_obj@meta.data)) {{
                p4 <- VlnPlot(seurat_obj, features = "percent.ribo", pt.size = 0.1) +
                      theme(legend.position = "none") +
                      labs(title = "Ribosomal Gene %", x = "", y = "Percent Ribosomal")
            }} else {{
                p4 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No ribosomal genes found") +
                      labs(title = "Ribosomal Gene %") +
                      theme_minimal()
            }}
            
            # Combine violin plots
            violin_combined <- p1 | p2 | p3 | p4
            
            ggsave("{output_dir}/qc_violin_plots.pdf", violin_combined, 
                   width = 16, height = 4, dpi = 300)
            
            cat("✅ Violin plots saved\\n")
            '''
            
            ro.r(violin_plot_code)
            
        except Exception as e:
            print(f"⚠️  Error generating violin plots: {e}")
        
        # Generate scatter plots with error handling
        try:
            scatter_plot_code = f'''
            # Scatter plot: UMI counts vs genes detected
            p1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
                  geom_smooth(method = "lm") +
                  labs(title = "UMI Counts vs Genes Detected")
            
            # Only create MT plots if MT data exists
            if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {{
                p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
                      geom_smooth(method = "lm") +
                      labs(title = "UMI Counts vs Mitochondrial %")
                
                p3 <- FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
                      geom_smooth(method = "lm") +
                      labs(title = "Genes Detected vs Mitochondrial %")
                
                scatter_combined <- p1 | p2 | p3
            }} else {{
                # Create placeholder plots if no MT data
                p2 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No MT data available") +
                      labs(title = "UMI Counts vs Mitochondrial %") +
                      theme_minimal()
                      
                p3 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No MT data available") +
                      labs(title = "Genes Detected vs Mitochondrial %") +
                      theme_minimal()
                      
                scatter_combined <- p1 | p2 | p3
            }}
            
            ggsave("{output_dir}/qc_scatter_plots.pdf", scatter_combined, 
                   width = 18, height = 6, dpi = 300)
            
            cat("✅ Scatter plots saved\\n")
            '''
            
            ro.r(scatter_plot_code)
            
        except Exception as e:
            print(f"⚠️  Error generating scatter plots: {e}")
        
        # Generate density plots with metadata directly
        try:
            density_plot_code = f'''
            # Get metadata for density plots
            metadata <- seurat_obj@meta.data
            
            # Density plot for UMI counts
            p1 <- ggplot(metadata, aes(x = nCount_RNA)) +
                  geom_density(fill = "steelblue", alpha = 0.7) +
                  scale_x_log10() +
                  labs(title = "UMI Counts Distribution", x = "UMI Counts (log10)", y = "Density") +
                  theme_minimal()
            
            # Density plot for genes detected
            p2 <- ggplot(metadata, aes(x = nFeature_RNA)) +
                  geom_density(fill = "darkgreen", alpha = 0.7) +
                  scale_x_log10() +
                  labs(title = "Genes Detected Distribution", x = "Genes Detected (log10)", y = "Density") +
                  theme_minimal()
            
            # Conditional MT and ribosomal plots
            if ("percent.mt" %in% colnames(metadata)) {{
                p3 <- ggplot(metadata, aes(x = percent.mt)) +
                      geom_density(fill = "salmon", alpha = 0.7) +
                      labs(title = "Mitochondrial % Distribution", x = "Percent MT", y = "Density") +
                      theme_minimal()
            }} else {{
                p3 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No MT data") +
                      labs(title = "Mitochondrial % Distribution") +
                      theme_minimal()
            }}
            
            if ("percent.ribo" %in% colnames(metadata)) {{
                p4 <- ggplot(metadata, aes(x = percent.ribo)) +
                      geom_density(fill = "orange", alpha = 0.7) +
                      labs(title = "Ribosomal % Distribution", x = "Percent Ribosomal", y = "Density") +
                      theme_minimal()
            }} else {{
                p4 <- ggplot() + 
                      annotate("text", x = 1, y = 1, label = "No ribosomal data") +
                      labs(title = "Ribosomal % Distribution") +
                      theme_minimal()
            }}
            
            density_combined <- (p1 | p2) / (p3 | p4)
            
            ggsave("{output_dir}/qc_density_plots.pdf", density_combined, 
                   width = 12, height = 8, dpi = 300)
            
            cat("✅ Density plots saved\\n")
            '''
            
            ro.r(density_plot_code)
            
        except Exception as e:
            print(f"⚠️  Error generating density plots: {e}")
        
        # Generate complexity plot if possible
        try:
            complexity_plot_code = f'''
            if ("log10GenesPerUMI" %in% colnames(metadata)) {{
                p_complexity <- ggplot(metadata, aes(x = log10GenesPerUMI)) +
                               geom_density(fill = "purple", alpha = 0.7) +
                               labs(title = "Complexity Score Distribution", 
                                    x = "Log10 Genes per UMI", y = "Density") +
                               theme_minimal()
                
                ggsave("{output_dir}/qc_complexity_plot.pdf", p_complexity, 
                       width = 8, height = 6, dpi = 300)
                
                cat("✅ Complexity plot saved\\n")
            }} else {{
                cat("⚠️  Skipping complexity plot - data not available\\n")
            }}
            '''
            
            ro.r(complexity_plot_code)
            
        except Exception as e:
            print(f"⚠️  Error generating complexity plot: {e}")
        
        print(f"✅ QC plots saved to: {output_dir}/")
        return True
    
    def detect_outliers_and_filter(self, seurat_obj_path: str, output_dir: str,
                                  min_genes: int = 200,
                                  max_genes: int = 5000,
                                  min_umis: int = 500,
                                  max_umis: int = 50000,
                                  max_mt_percent: float = 20.0,
                                  min_complexity: float = 0.8):
        """
        Detect outliers and apply filtering based on tutorial recommendations
        """
        print("🎯 Detecting outliers and filtering cells...")
        
        # Load Seurat object
        ro.r(f'seurat_obj <- readRDS("{seurat_obj_path}")')
        
        # Get pre-filter statistics - check if QC metrics exist first
        pre_filter_code = '''
        # Check what columns are available
        available_cols <- colnames(seurat_obj@meta.data)
        cat("Available columns in metadata:\\n")
        cat(paste(available_cols, collapse = ", "), "\\n\\n")
        
        # Get basic statistics
        pre_filter_stats <- list(
            n_cells = ncol(seurat_obj),
            n_genes = nrow(seurat_obj),
            median_genes_per_cell = median(seurat_obj$nFeature_RNA),
            median_umis_per_cell = median(seurat_obj$nCount_RNA)
        )
        
        # Add MT stats if available
        if ("percent.mt" %in% available_cols) {
            pre_filter_stats$median_mt_percent <- median(seurat_obj$percent.mt)
            cat("✅ Mitochondrial data available\\n")
        } else {
            cat("❌ No mitochondrial data found\\n")
        }
        
        # Add complexity stats if available  
        if ("log10GenesPerUMI" %in% available_cols) {
            pre_filter_stats$median_complexity <- median(seurat_obj$log10GenesPerUMI, na.rm = TRUE)
            cat("✅ Complexity data available\\n")
        } else {
            cat("❌ No complexity data found\\n")
        }
        '''
        
        ro.r(pre_filter_code)
        
        # Apply filters step by step with reporting
        filter_code = f'''
        # Ensure QC metrics exist before filtering
        if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {{
            seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
        }}
        if (!"log10GenesPerUMI" %in% colnames(seurat_obj@meta.data)) {{
            seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
        }}
        
        # Initial cell count
        cat("Initial number of cells:", ncol(seurat_obj), "\\n")
        
        # Filter 1: Remove cells with too few or too many genes
        cells_keep_genes <- seurat_obj$nFeature_RNA >= {min_genes} & seurat_obj$nFeature_RNA <= {max_genes}
        cat("Cells passing gene count filter:", sum(cells_keep_genes), "\\n")
        cat("Cells removed by gene count filter:", sum(!cells_keep_genes), "\\n")
        
        # Filter 2: Remove cells with too few or too many UMIs
        cells_keep_umis <- seurat_obj$nCount_RNA >= {min_umis} & seurat_obj$nCount_RNA <= {max_umis}
        cat("Cells passing UMI count filter:", sum(cells_keep_umis), "\\n")
        cat("Cells removed by UMI count filter:", sum(!cells_keep_umis), "\\n")
        
        # Filter 3: Remove cells with high mitochondrial content (if data available)
        if (mean(seurat_obj$percent.mt, na.rm = TRUE) > 0) {{
            cells_keep_mt <- seurat_obj$percent.mt <= {max_mt_percent}
            cat("Cells passing MT% filter:", sum(cells_keep_mt), "\\n")
            cat("Cells removed by MT% filter:", sum(!cells_keep_mt), "\\n")
        }} else {{
            cells_keep_mt <- rep(TRUE, ncol(seurat_obj))
            cat("No mitochondrial data available - skipping MT filter\\n")
        }}
        
        # Filter 4: Remove cells with low complexity (if data valid)
        if (all(is.finite(seurat_obj$log10GenesPerUMI))) {{
            cells_keep_complexity <- seurat_obj$log10GenesPerUMI >= {min_complexity}
            cat("Cells passing complexity filter:", sum(cells_keep_complexity), "\\n")
            cat("Cells removed by complexity filter:", sum(!cells_keep_complexity), "\\n")
        }} else {{
            cells_keep_complexity <- rep(TRUE, ncol(seurat_obj))
            cat("Invalid complexity scores - skipping complexity filter\\n")
        }}
        
        # Combine all filters
        cells_to_keep <- cells_keep_genes & cells_keep_umis & cells_keep_mt & cells_keep_complexity
        
        cat("\\nFinal filtering results:\\n")
        cat("Cells before filtering:", ncol(seurat_obj), "\\n")
        cat("Cells after filtering:", sum(cells_to_keep), "\\n")
        cat("Cells removed:", sum(!cells_to_keep), "\\n")
        cat("Percentage of cells kept:", round(sum(cells_to_keep)/ncol(seurat_obj)*100, 2), "%\\n")
        
        # Apply filter
        seurat_filtered <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
        
        # Gene filtering: Remove genes expressed in fewer than 3 cells
        genes_to_keep <- rowSums(GetAssayData(seurat_filtered, slot = "counts") > 0) >= 3
        seurat_filtered <- seurat_filtered[genes_to_keep, ]
        
        cat("\\nGene filtering results:\\n")
        cat("Genes before filtering:", nrow(seurat_obj), "\\n")
        cat("Genes after filtering:", nrow(seurat_filtered), "\\n")
        cat("Genes removed:", nrow(seurat_obj) - nrow(seurat_filtered), "\\n")
        '''
        
        ro.r(filter_code)
        
        # Save filtered object
        filtered_path = os.path.join(output_dir, "seurat_filtered.rds")
        ro.r(f'saveRDS(seurat_filtered, "{filtered_path}")')
        
        # Generate comparison plots with error handling
        try:
            comparison_plot_code = f'''
            # Create before/after comparison plots
            
            # Prepare data for comparison
            metadata_before <- seurat_obj@meta.data
            metadata_before$Filter_Status <- "Before"
            metadata_after <- seurat_filtered@meta.data
            metadata_after$Filter_Status <- "After"
            
            # Select common columns that exist in both datasets
            common_cols <- c("nFeature_RNA", "nCount_RNA", "Filter_Status")
            if ("percent.mt" %in% colnames(metadata_before) && "percent.mt" %in% colnames(metadata_after)) {{
                common_cols <- c(common_cols, "percent.mt")
            }}
            
            # Combine metadata
            combined_metadata <- rbind(
                metadata_before[, common_cols],
                metadata_after[, common_cols]
            )
            
            # Comparison violin plots
            p1 <- ggplot(combined_metadata, aes(x = Filter_Status, y = nFeature_RNA, fill = Filter_Status)) +
                  geom_violin(alpha = 0.7) +
                  geom_boxplot(width = 0.1, fill = "white") +
                  labs(title = "Genes Detected: Before vs After Filtering", 
                       x = "", y = "Number of Genes") +
                  theme_minimal() +
                  theme(legend.position = "none")
            
            p2 <- ggplot(combined_metadata, aes(x = Filter_Status, y = nCount_RNA, fill = Filter_Status)) +
                  geom_violin(alpha = 0.7) +
                  geom_boxplot(width = 0.1, fill = "white") +
                  labs(title = "UMI Counts: Before vs After Filtering", 
                       x = "", y = "UMI Counts") +
                  theme_minimal() +
                  theme(legend.position = "none")
            
            # Only add MT plot if data exists
            if ("percent.mt" %in% common_cols) {{
                p3 <- ggplot(combined_metadata, aes(x = Filter_Status, y = percent.mt, fill = Filter_Status)) +
                      geom_violin(alpha = 0.7) +
                      geom_boxplot(width = 0.1, fill = "white") +
                      labs(title = "MT%: Before vs After Filtering", 
                           x = "", y = "Percent MT") +
                      theme_minimal() +
                      theme(legend.position = "none")
                
                comparison_plots <- p1 | p2 | p3
            }} else {{
                comparison_plots <- p1 | p2
            }}
            
            ggsave("{output_dir}/filtering_comparison.pdf", comparison_plots, 
                   width = 15, height = 5, dpi = 300)
            
            cat("✅ Comparison plots saved\\n")
            '''
            
            ro.r(comparison_plot_code)
            
        except Exception as e:
            print(f"⚠️  Error generating comparison plots: {e}")
        
        print(f"✅ Filtered Seurat object saved to: {filtered_path}")
        return filtered_path
    
    def generate_comprehensive_qc_report(self, h5ad_path: str, output_dir: str = "./comprehensive_qc_output"):
        """
        Generate a complete QC report following the tutorial workflow
        """
        print("=" * 60)
        print("🧬 COMPREHENSIVE SINGLE CELL QC ANALYSIS")
        print("=" * 60)
        print(f"Input H5AD file: {h5ad_path}")
        print(f"Output directory: {output_dir}")
        print("")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Step 1: Convert H5AD to Seurat RDS
        rds_path = os.path.join(output_dir, "seurat_object.rds")
        seurat_obj, adata = self.convert_h5ad_to_rds(h5ad_path, rds_path, 
                                                    save_coords=True, normalize=False)
        
        # Step 2: Calculate QC metrics
        metadata = self.calculate_qc_metrics(rds_path)
        
        # Step 3: Print summary statistics
        print("\\n📈 QC METRICS SUMMARY:")
        print("-" * 40)
        qc_metrics = ['nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo']
        for metric in qc_metrics:
            if metric in metadata.columns:
                values = metadata[metric]
                print(f"{metric}:")
                print(f"  Mean ± SD: {values.mean():.1f} ± {values.std():.1f}")
                print(f"  Median (IQR): {values.median():.1f} ({values.quantile(0.25):.1f}-{values.quantile(0.75):.1f})")
                print(f"  Range: {values.min():.1f} - {values.max():.1f}")
                print("")
        
        # Step 4: Generate QC plots
        self.generate_qc_plots_tutorial_style(rds_path, output_dir)
        
        # Step 5: Detect outliers and filter cells
        filtered_rds_path = self.detect_outliers_and_filter(rds_path, output_dir)
        
        # Step 6: Generate final summary
        print("\\n📊 FINAL SUMMARY:")
        print("-" * 30)
        ro.r(f'original_obj <- readRDS("{rds_path}")')
        ro.r(f'filtered_obj <- readRDS("{filtered_rds_path}")')
        
        summary_code = '''
        cat("Original dataset:\\n")
        cat("  Cells:", ncol(original_obj), "\\n")
        cat("  Genes:", nrow(original_obj), "\\n")
        cat("\\nFiltered dataset:\\n")
        cat("  Cells:", ncol(filtered_obj), "\\n")
        cat("  Genes:", nrow(filtered_obj), "\\n")
        cat("\\nFiltering efficiency:\\n")
        cat("  Cells retained:", round(ncol(filtered_obj)/ncol(original_obj)*100, 2), "%\\n")
        cat("  Genes retained:", round(nrow(filtered_obj)/nrow(original_obj)*100, 2), "%\\n")
        '''
        ro.r(summary_code)
        
        # Step 7: Save metadata
        metadata_path = os.path.join(output_dir, "qc_metadata.csv")
        metadata.to_csv(metadata_path)
        print(f"\\n💾 Metadata saved to: {metadata_path}")
        
        print("\\n✅ COMPREHENSIVE QC ANALYSIS COMPLETE!")
        print(f"📁 All results saved to: {output_dir}")
        
        return {
            'original_rds': rds_path,
            'filtered_rds': filtered_rds_path,
            'metadata': metadata,
            'output_dir': output_dir,
            'adata': adata
        }

    def generate_comprehensive_qc_report_v2(self, h5ad_path: str, output_dir: str = "./comprehensive_qc_output",
                                            # QC filtering parameters following PBMC3k tutorial
                                            min_features: int = 200,
                                            max_features: int = 4000, 
                                            max_mt_percent: float = 5.0,
                                            # Advanced filtering options
                                            use_advanced_filtering: bool = True,
                                            min_umis: int = 500,
                                            max_umis: int = 50000,
                                            min_complexity: float = 0.8,
                                            # Analysis parameters
                                            n_variable_features: int = 2000,
                                            n_pcs: int = 10,
                                            clustering_resolution: float = 0.5,
                                            # Advanced options
                                            regress_out_vars: List[str] = None,
                                            use_sct_transform: bool = False):
        """
        Generate a complete QC report following the official Seurat PBMC3k tutorial workflow
        
        This function implements the standard Seurat preprocessing pipeline:
        1. Data loading and QC metric calculation
        2. Quality control filtering (following PBMC3k standards)
        3. Data normalization (LogNormalize or SCTransform)
        4. Highly variable feature detection
        5. Data scaling and PCA
        6. Clustering and UMAP generation
        7. Marker gene identification
        8. Comprehensive visualization and reporting
        
        Parameters:
        -----------
        h5ad_path : str
            Path to input H5AD file
        output_dir : str
            Output directory for results
        min_features : int
            Minimum number of features (genes) per cell (default: 200, as in tutorial)
        max_features : int
            Maximum number of features per cell (default: 2500, as in tutorial)
        max_mt_percent : float
            Maximum mitochondrial gene percentage (default: 5.0, as in tutorial)
        n_variable_features : int
            Number of highly variable features to find (default: 2000, as in tutorial)
        n_pcs : int
            Number of principal components to use (default: 10, as in tutorial)
        clustering_resolution : float
            Clustering resolution (default: 0.5, as in tutorial)
        regress_out_vars : List[str]
            Variables to regress out during scaling (e.g., ['percent.mt'])
        use_sct_transform : bool
            Whether to use SCTransform instead of standard normalization
        
        Returns:
        --------
        dict: Dictionary containing analysis results and file paths
        """
        print("=" * 80)
        print("🧬 COMPREHENSIVE SINGLE CELL QC ANALYSIS")
        print("📚 Following Seurat PBMC3k Tutorial Pipeline")
        print("🔗 https://satijalab.org/seurat/articles/pbmc3k_tutorial.html")
        print("=" * 80)
        print(f"📁 Input H5AD file: {h5ad_path}")
        print(f"📁 Output directory: {output_dir}")
        print(f"🎯 QC Parameters: min_features={min_features}, max_features={max_features}, max_mt={max_mt_percent}%")
        print("")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # =====================================
        # STEP 1: DATA LOADING AND CONVERSION
        # =====================================
        print("📊 STEP 1: DATA LOADING AND SEURAT OBJECT CREATION")
        print("-" * 60)
        
        rds_path = os.path.join(output_dir, "seurat_object_raw.rds")
        seurat_obj, adata = self.convert_h5ad_to_rds(h5ad_path, rds_path, 
                                                    save_coords=True, normalize=False)
        
        # =====================================
        # STEP 2: QC METRICS CALCULATION
        # =====================================
        print("\n🔬 STEP 2: QC METRICS CALCULATION")
        print("-" * 60)
        
        # Calculate comprehensive QC metrics following tutorial
        ro.r(f'seurat_obj <- readRDS("{rds_path}")')
        
        qc_metrics_code = '''
        library(Seurat)
        
        # Calculate mitochondrial gene percentage (main QC metric in tutorial)
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Mt-")
        
        # Additional QC metrics for comprehensive analysis
        seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
        seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(p)]")
        seurat_obj[["percent.platelet"]] <- PercentageFeatureSet(seurat_obj, pattern = "PPBP|PF4|SDPR|GP1BA|TUBB1|CLU|GNG11|NRGN|RGS18|TPM4")
        
        # Calculate complexity score (genes per UMI)
        seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
        
        # Print QC summary as in tutorial
        cat("📈 QC METRICS SUMMARY (Raw Data):\\n")
        cat("================================\\n")
        cat("Number of cells:", ncol(seurat_obj), "\\n")
        cat("Number of genes:", nrow(seurat_obj), "\\n")
        cat("Median genes per cell:", median(seurat_obj$nFeature_RNA), "\\n")
        cat("Median UMI per cell:", median(seurat_obj$nCount_RNA), "\\n")
        cat("Median MT% per cell:", round(median(seurat_obj$percent.mt), 2), "%\\n")
        cat("\\n")
        '''
        
        ro.r(qc_metrics_code)
        
        # Save object with QC metrics
        ro.r(f'saveRDS(seurat_obj, "{rds_path}")')
        
        # =====================================
        # STEP 3: QC VISUALIZATION (Pre-filtering)
        # =====================================
        print("📊 STEP 3: QC VISUALIZATION (Pre-filtering)")
        print("-" * 60)
        
        # Generate pre-filtering QC plots following tutorial style
        pre_filter_plots_code = f'''
        # Create QC plots as shown in PBMC3k tutorial
        library(patchwork)
        
        # Violin plots for main QC metrics (exactly as in tutorial)
        p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
                    pt.size = 0.1) & theme(plot.title = element_text(size=10))
        
        ggsave("{output_dir}/01_qc_violin_plots_raw.pdf", p1, width = 12, height = 4, dpi = 300)
        
        # Feature scatter plots (as in tutorial)
        plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
                geom_smooth(method = "lm", se = FALSE)
        plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
                geom_smooth(method = "lm", se = FALSE)
        
        combined_scatter <- plot1 + plot2
        ggsave("{output_dir}/02_qc_scatter_plots_raw.pdf", combined_scatter, width = 12, height = 5, dpi = 300)
        
        # Additional density plots for better understanding
        metadata <- seurat_obj@meta.data
        
        p_density <- ggplot(metadata, aes(x = nFeature_RNA)) +
                    geom_density(fill = "steelblue", alpha = 0.7) +
                    geom_vline(xintercept = {min_features}, color = "red", linetype = "dashed") +
                    geom_vline(xintercept = {max_features}, color = "red", linetype = "dashed") +
                    labs(title = "Gene Count Distribution", x = "Number of Genes", y = "Density") +
                    theme_minimal()
        
        p_mt <- ggplot(metadata, aes(x = percent.mt)) +
            geom_density(fill = "salmon", alpha = 0.7) +
            geom_vline(xintercept = {max_mt_percent}, color = "red", linetype = "dashed") +
            labs(title = "Mitochondrial % Distribution", x = "Percent MT", y = "Density") +
            theme_minimal()
        
        p_combined <- p_density | p_mt
        ggsave("{output_dir}/03_qc_density_plots_raw.pdf", p_combined, width = 12, height = 4, dpi = 300)
        
        cat("✅ Pre-filtering QC plots saved\\n")
        '''
        
        ro.r(pre_filter_plots_code)
        
        # =====================================
        # STEP 4: CELL FILTERING
        # =====================================
        print("\n🎯 STEP 4: CELL AND GENE FILTERING")
        print("-" * 60)
        
        # Option 1: Use the comprehensive outlier detection function
        if use_advanced_filtering:
            print("🔬 Using comprehensive outlier detection and filtering...")
            filtered_path = self.detect_outliers_and_filter(
                rds_path, output_dir,
                min_genes=min_features,
                max_genes=max_features,
                min_umis=min_umis,
                max_umis=max_umis,
                max_mt_percent=max_mt_percent,
                min_complexity=min_complexity
            )
            ro.r(f'seurat_filtered <- readRDS("{filtered_path}")')
            
        else:
            # Option 2: Simple filtering following PBMC3k tutorial exactly
            print("📊 Using PBMC3k tutorial standard filtering...")
            filtering_code = f'''
            # Cell filtering following PBMC3k tutorial exactly
            cat("Filtering cells based on PBMC3k tutorial criteria:\\n")
            cat("- nFeature_RNA > {min_features} & nFeature_RNA < {max_features}\\n")
            cat("- percent.mt < {max_mt_percent}\\n\\n")
            
            # Show before filtering stats
            cat("Before filtering:\\n")
            cat("  Cells:", ncol(seurat_obj), "\\n")
            cat("  Genes:", nrow(seurat_obj), "\\n")
            
            # Apply the same filtering as PBMC3k tutorial
            seurat_filtered <- subset(seurat_obj, subset = nFeature_RNA > {min_features} & 
                                                        nFeature_RNA < {max_features} & 
                                                        percent.mt < {max_mt_percent})
            
            # Show after filtering stats
            cat("\\nAfter filtering:\\n")
            cat("  Cells:", ncol(seurat_filtered), "\\n")
            cat("  Genes:", nrow(seurat_filtered), "\\n")
            cat("  Cells retained:", round(ncol(seurat_filtered)/ncol(seurat_obj)*100, 2), "%\\n\\n")
            '''
            
            ro.r(filtering_code)
            
            # Save filtered object
            filtered_path = os.path.join(output_dir, "seurat_object_filtered.rds")
            ro.r(f'saveRDS(seurat_filtered, "{filtered_path}")')
            
        # =====================================
        # STEP 5: DATA NORMALIZATION
        # =====================================
        print("\n🔄 STEP 5: DATA NORMALIZATION")
        print("-" * 60)
        
        if use_sct_transform:
            print("📊 Using SCTransform normalization...")
            normalize_code = f'''
            # SCTransform normalization (alternative method)
            seurat_filtered <- SCTransform(seurat_filtered, verbose = FALSE)
            cat("✅ SCTransform normalization completed\\n")
            '''
        else:
            print("📊 Using standard LogNormalize (as in PBMC3k tutorial)...")
            normalize_code = f'''
            # Standard normalization as in PBMC3k tutorial
            seurat_filtered <- NormalizeData(seurat_filtered, normalization.method = "LogNormalize", 
                                        scale.factor = 10000, verbose = FALSE)
            cat("✅ LogNormalize normalization completed\\n")
            '''
        
        ro.r(normalize_code)
        
        # =====================================
        # STEP 6: HIGHLY VARIABLE FEATURES
        # =====================================
        print("\n🔍 STEP 6: IDENTIFICATION OF HIGHLY VARIABLE FEATURES")
        print("-" * 60)
        
        if not use_sct_transform:  # SCTransform includes variable feature detection
            hvf_code = f'''
            # Find variable features as in PBMC3k tutorial
            seurat_filtered <- FindVariableFeatures(seurat_filtered, selection.method = "vst", 
                                                nfeatures = {n_variable_features}, verbose = FALSE)
            
            # Identify the 10 most highly variable genes
            top10 <- head(VariableFeatures(seurat_filtered), 10)
            
            # Plot variable features with and without labels (as in tutorial)
            plot1 <- VariableFeaturePlot(seurat_filtered)
            plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
            combined_hvf <- plot1 + plot2
            
            ggsave("{output_dir}/04_highly_variable_features.pdf", combined_hvf, 
                width = 15, height = 6, dpi = 300)
            
            cat("✅ Highly variable features identified:", length(VariableFeatures(seurat_filtered)), "\\n")
            cat("Top 10 HVFs:", paste(top10, collapse = ", "), "\\n")
            '''
            ro.r(hvf_code)
        
        # =====================================
        # STEP 7: DATA SCALING AND PCA
        # =====================================
        print("\n📐 STEP 7: DATA SCALING AND PRINCIPAL COMPONENT ANALYSIS")
        print("-" * 60)
        
        if not use_sct_transform:
            if regress_out_vars:
                print(f"🔧 Regressing out variables: {regress_out_vars}")
                scale_vars = ', '.join([f'"{var}"' for var in regress_out_vars])
                scaling_code = f'''
                # Scale data with regression (optional, as mentioned in tutorial)
                seurat_filtered <- ScaleData(seurat_filtered, vars.to.regress = c({scale_vars}), verbose = FALSE)
                '''
            else:
                scaling_code = f'''
                # Standard scaling as in tutorial
                seurat_filtered <- ScaleData(seurat_filtered, verbose = FALSE)
                '''
            
            ro.r(scaling_code)
        
        # PCA analysis
        pca_code = f'''
        # PCA as in PBMC3k tutorial
        seurat_filtered <- RunPCA(seurat_filtered, features = VariableFeatures(object = seurat_filtered), 
                                verbose = FALSE)
        
        # Print PCA results as in tutorial
        print(seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
        
        # Visualize PCA results
        p1 <- VizDimLoadings(seurat_filtered, dims = 1:2, reduction = "pca", ncol = 2)
        ggsave("{output_dir}/05_pca_loadings.pdf", p1, width = 10, height = 8, dpi = 300)
        
        # PCA plot
        p2 <- DimPlot(seurat_filtered, reduction = "pca") + ggtitle("PCA")
        ggsave("{output_dir}/06_pca_plot.pdf", p2, width = 8, height = 6, dpi = 300)
        
        # Elbow plot for PC selection (as in tutorial)
        p3 <- ElbowPlot(seurat_filtered, ndims = 20)
        ggsave("{output_dir}/07_elbow_plot.pdf", p3, width = 8, height = 6, dpi = 300)
        
        # Heatmap of PCs (as in tutorial)
        p4 <- DimHeatmap(seurat_filtered, dims = 1:6, cells = 500, balanced = TRUE, ncol = 2)
        ggsave("{output_dir}/08_pca_heatmap.pdf", p4, width = 12, height = 12, dpi = 300)
        
        cat("✅ PCA analysis completed\\n")
        '''
        
        ro.r(pca_code)
        
        # =====================================
        # STEP 8: CLUSTERING
        # =====================================
        print("\n🔗 STEP 8: CLUSTERING ANALYSIS")
        print("-" * 60)
        
        clustering_code = f'''
        # Clustering following PBMC3k tutorial exactly
        seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:{n_pcs}, verbose = FALSE)
        seurat_filtered <- FindClusters(seurat_filtered, resolution = {clustering_resolution}, verbose = FALSE)
        
        # Check cluster assignments
        cat("Number of clusters found:", length(unique(Idents(seurat_filtered))), "\\n")
        cluster_table <- table(Idents(seurat_filtered))
        cat("Cells per cluster:\\n")
        print(cluster_table)
        '''
        
        ro.r(clustering_code)
        
        # =====================================
        # STEP 9: NON-LINEAR DIMENSIONALITY REDUCTION
        # =====================================
        print("\n🎨 STEP 9: UMAP VISUALIZATION")
        print("-" * 60)
        
        umap_code = f'''
        # UMAP as in PBMC3k tutorial
        seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:{n_pcs}, verbose = FALSE)
        
        # Create UMAP plots
        p1 <- DimPlot(seurat_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + 
            ggtitle("UMAP - Clusters") + NoLegend()
        
        # Also create unlabeled version
        p2 <- DimPlot(seurat_filtered, reduction = "umap", pt.size = 0.5) + 
            ggtitle("UMAP - Clusters (with legend)")
        
        combined_umap <- p1 | p2
        ggsave("{output_dir}/09_umap_clusters.pdf", combined_umap, width = 15, height = 6, dpi = 300)
        
        # QC metrics on UMAP
        p3 <- FeaturePlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        reduction = "umap", ncol = 3)
        ggsave("{output_dir}/10_umap_qc_features.pdf", p3, width = 15, height = 5, dpi = 300)
        
        cat("✅ UMAP visualization completed\\n")
        '''
        
        ro.r(umap_code)
        
        # =====================================
        # STEP 10: MARKER GENE IDENTIFICATION
        # =====================================
        print("\n🧬 STEP 10: MARKER GENE IDENTIFICATION")
        print("-" * 60)
        
        markers_code = f'''
        # Find markers for all clusters (as in tutorial)
        all_markers <- FindAllMarkers(seurat_filtered, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25, verbose = FALSE)
        
        # Get top markers per cluster
        top_markers <- all_markers %>% 
                    group_by(cluster) %>% 
                    filter(avg_log2FC > 1) %>%
                    slice_head(n = 10)
        
        cat("✅ Marker genes identified for", length(unique(all_markers$cluster)), "clusters\\n")
        cat("Total significant markers:", nrow(all_markers), "\\n")
        
        # Save marker results
        write.csv(all_markers, "{output_dir}/11_all_cluster_markers.csv", row.names = FALSE)
        write.csv(top_markers, "{output_dir}/12_top_cluster_markers.csv", row.names = FALSE)
        '''
        
        ro.r(markers_code)
        
        # =====================================
        # STEP 11: VISUALIZATION OF MARKER GENES
        # =====================================
        print("\n📊 STEP 11: MARKER GENE VISUALIZATION")
        print("-" * 60)
        
        viz_code = f'''
        # Create marker visualization plots
        
        # Heatmap of top markers
        top10_markers <- all_markers %>% 
                        group_by(cluster) %>% 
                        filter(avg_log2FC > 1) %>%
                        slice_head(n = 10) %>% 
                        ungroup()
        
        if(nrow(top10_markers) > 0) {{
            p1 <- DoHeatmap(seurat_filtered, features = top10_markers$gene) + NoLegend()
            ggsave("{output_dir}/13_marker_heatmap.pdf", p1, width = 12, height = 15, dpi = 300)
        }}
        
        # Get some example markers for feature plots (if they exist)
        common_markers <- c("Syt1", "Plp1", "Cldn5", "Cq1a", "Pdgfra", "Gad1", "Slc17a7", "Aqp4")
        available_markers <- intersect(common_markers, rownames(seurat_filtered))
        
        if(length(available_markers) > 0) {{
            p2 <- FeaturePlot(seurat_filtered, features = available_markers[1:min(9, length(available_markers))], 
                            ncol = 3, pt.size = 0.3)
            ggsave("{output_dir}/14_common_marker_features.pdf", p2, width = 15, height = 15, dpi = 300)
            
            # Violin plots for some markers
            p3 <- VlnPlot(seurat_filtered, features = available_markers[1:min(4, length(available_markers))], 
                        ncol = 2, pt.size = 0.1)
            ggsave("{output_dir}/15_marker_violin_plots.pdf", p3, width = 12, height = 8, dpi = 300)
        }}
        
        cat("✅ Marker gene visualizations created\\n")
        '''
        
        ro.r(viz_code)
        
        # =====================================
        # STEP 12: GENERATE COMPREHENSIVE REPORT
        # =====================================
        print("\n📋 STEP 12: GENERATING FINAL SUMMARY REPORT")
        print("-" * 60)
        
        # # Generate final summary and comparison
        # summary_code = f'''
        # # Final summary statistics
        # cat("\\n" %R% "=" * 60 %R% "\\n")
        # cat("📊 FINAL ANALYSIS SUMMARY\\n")
        # cat("=" * 60 %R% "\\n")
        
        # # Original vs filtered comparison
        # original_obj <- readRDS("{rds_path}")
        
        # cat("Original dataset (after initial processing):\\n")
        # cat("  Cells:", ncol(original_obj), "\\n")
        # cat("  Genes:", nrow(original_obj), "\\n")
        # cat("  Median genes/cell:", median(original_obj$nFeature_RNA), "\\n")
        # cat("  Median UMI/cell:", median(original_obj$nCount_RNA), "\\n")
        # cat("  Median MT%:", round(median(original_obj$percent.mt), 2), "%\\n\\n")
        
        # cat("Filtered and processed dataset:\\n")
        # cat("  Cells:", ncol(seurat_filtered), "\\n")
        # cat("  Genes:", nrow(seurat_filtered), "\\n")
        # cat("  Median genes/cell:", median(seurat_filtered$nFeature_RNA), "\\n")
        # cat("  Median UMI/cell:", median(seurat_filtered$nCount_RNA), "\\n")
        # cat("  Median MT%:", round(median(seurat_filtered$percent.mt), 2), "%\\n")
        # cat("  Number of clusters:", length(unique(Idents(seurat_filtered))), "\\n")
        # cat("  Number of HVGs:", length(VariableFeatures(seurat_filtered)), "\\n\\n")
        
        # cat("Filtering efficiency:\\n")
        # cat("  Cells retained:", round(ncol(seurat_filtered)/ncol(original_obj)*100, 2), "%\\n")
        # cat("  Genes retained:", round(nrow(seurat_filtered)/nrow(original_obj)*100, 2), "%\\n\\n")
        
        # # Analysis parameters used
        # cat("Analysis parameters used:\\n")
        # cat("  Min features per cell:", {min_features}, "\\n")
        # cat("  Max features per cell:", {max_features}, "\\n") 
        # cat("  Max mitochondrial %:", {max_mt_percent}, "%\\n")
        # cat("  Normalization method:", "{"SCTransform" if use_sct_transform else "LogNormalize"}", "\\n")
        # cat("  Variable features:", {n_variable_features}, "\\n")
        # cat("  PCs used for clustering:", {n_pcs}, "\\n")
        # cat("  Clustering resolution:", {clustering_resolution}, "\\n")
        # '''
        
        # # Fix R string concatenation
        # summary_code = summary_code.replace('%R%', '')
        # summary_code = summary_code.replace('"=" * 60', 'paste(rep("=", 60), collapse="")')
        
        # ro.r(summary_code)
        # =====================================
        # STEP 12: GENERATE COMPREHENSIVE REPORT
        # =====================================
        print("\n📋 STEP 12: GENERATING FINAL SUMMARY REPORT")
        print("-" * 60)
        
        # Generate final summary and comparison
        print("Generating final summary statistics...")
        
        # Load original object for comparison
        ro.r(f'original_obj <- readRDS("{rds_path}")')
        
        # Generate summary with proper R syntax
        summary_code = '''
        # Final summary statistics
        separator_line <- paste(rep("=", 60), collapse="")
        cat("\\n", separator_line, "\\n", sep="")
        cat("FINAL ANALYSIS SUMMARY\\n")
        cat(separator_line, "\\n", sep="")
        
        cat("Original dataset (after initial processing):\\n")
        cat("  Cells:", ncol(original_obj), "\\n")
        cat("  Genes:", nrow(original_obj), "\\n")
        cat("  Median genes/cell:", median(original_obj$nFeature_RNA), "\\n")
        cat("  Median UMI/cell:", median(original_obj$nCount_RNA), "\\n")
        cat("  Median MT%:", round(median(original_obj$percent.mt), 2), "%\\n\\n")
        
        cat("Filtered and processed dataset:\\n")
        cat("  Cells:", ncol(seurat_filtered), "\\n")
        cat("  Genes:", nrow(seurat_filtered), "\\n")
        cat("  Median genes/cell:", median(seurat_filtered$nFeature_RNA), "\\n")
        cat("  Median UMI/cell:", median(seurat_filtered$nCount_RNA), "\\n")
        cat("  Median MT%:", round(median(seurat_filtered$percent.mt), 2), "%\\n")
        cat("  Number of clusters:", length(unique(Idents(seurat_filtered))), "\\n")
        cat("  Number of HVGs:", length(VariableFeatures(seurat_filtered)), "\\n\\n")
        
        cat("Filtering efficiency:\\n")
        cat("  Cells retained:", round(ncol(seurat_filtered)/ncol(original_obj)*100, 2), "%\\n")
        cat("  Genes retained:", round(nrow(seurat_filtered)/nrow(original_obj)*100, 2), "%\\n\\n")
        '''
        
        ro.r(summary_code)
    
        # Save final processed object
        final_path = os.path.join(output_dir, "seurat_object_final.rds")
        ro.r(f'saveRDS(seurat_filtered, "{final_path}")')
        
        # Get final metadata for return
        metadata_r = ro.r('seurat_filtered@meta.data')
        with localconverter(default_converter + pandas2ri.converter):
            final_metadata = pandas2ri.rpy2py(metadata_r)
        
        # Save metadata
        metadata_path = os.path.join(output_dir, "final_metadata.csv")
        final_metadata.to_csv(metadata_path)
        
        # =====================================
        # STEP 13: CREATE FILE SUMMARY
        # =====================================
        print("\n📁 STEP 13: FILE SUMMARY")
        print("-" * 60)
        
        # Create a summary of all generated files
        file_summary = f"""
                        # Single Cell RNA-seq Analysis Results
                        ## Following Seurat PBMC3k Tutorial Pipeline

                        ### Analysis Overview
                        - Input file: {h5ad_path}
                        - Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
                        - Pipeline: Seurat PBMC3k Tutorial Standard Workflow

                        ### Quality Control Parameters
                        - Minimum features per cell: {min_features}
                        - Maximum features per cell: {max_features}  
                        - Maximum mitochondrial percentage: {max_mt_percent}%

                        ### Analysis Parameters
                        - Normalization method: {"SCTransform" if use_sct_transform else "LogNormalize"}
                        - Variable features identified: {n_variable_features}
                        - Principal components used: {n_pcs}
                        - Clustering resolution: {clustering_resolution}
                        - Variables regressed out: {regress_out_vars if regress_out_vars else "None"}

                        ### Generated Files

                        #### Seurat Objects
                        - `seurat_object_raw.rds` - Initial object with QC metrics
                        - `seurat_object_filtered.rds` - After quality filtering
                        - `seurat_object_final.rds` - Final processed object with clustering

                        #### Quality Control Plots  
                        - `01_qc_violin_plots_raw.pdf` - Pre-filtering QC metrics
                        - `02_qc_scatter_plots_raw.pdf` - Feature-feature relationships
                        - `03_qc_density_plots_raw.pdf` - Distribution plots with filter thresholds

                        #### Analysis Plots
                        - `04_highly_variable_features.pdf` - Variable feature identification
                        - `05_pca_loadings.pdf` - PCA gene loadings
                        - `06_pca_plot.pdf` - PCA scatter plot
                        - `07_elbow_plot.pdf` - PC selection elbow plot  
                        - `08_pca_heatmap.pdf` - Heatmap of top PCs

                        #### Clustering and Visualization
                        - `09_umap_clusters.pdf` - UMAP with cluster labels
                        - `10_umap_qc_features.pdf` - QC metrics on UMAP

                        #### Marker Analysis
                        - `11_all_cluster_markers.csv` - Complete marker gene results
                        - `12_top_cluster_markers.csv` - Top 10 markers per cluster
                        - `13_marker_heatmap.pdf` - Heatmap of top markers
                        - `14_common_marker_features.pdf` - Feature plots of known markers
                        - `15_marker_violin_plots.pdf` - Violin plots of marker expression

                        #### Metadata
                        - `final_metadata.csv` - Cell metadata with cluster assignments
                        - `analysis_summary.md` - This summary file

                        ### Next Steps
                        1. Review QC plots to validate filtering parameters
                        2. Examine cluster markers for biological interpretation  
                        3. Consider cell type annotation using known markers
                        4. Perform downstream analyses (differential expression, pathway analysis)

                        ### Tutorial Reference
                        This analysis follows the official Seurat PBMC3k tutorial:
                        https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
                        """
        
        summary_path = os.path.join(output_dir, "analysis_summary.md")
        with open(summary_path, 'w') as f:
            f.write(file_summary)
        
        print("\n🎉 ANALYSIS COMPLETED SUCCESSFULLY!")
        print("=" * 80)
        print(f"📁 All results saved to: {output_dir}")
        print(f"📋 Analysis summary: {summary_path}")
        print(f"🧬 Final Seurat object: {final_path}")
        print(f"📊 Metadata file: {metadata_path}")
        
        # Return comprehensive results
        results = {
            'input_file': h5ad_path,
            'output_directory': output_dir,
            'seurat_objects': {
                'raw': rds_path,
                'filtered': filtered_path,
                'final': final_path
            },
            'metadata_file': metadata_path,
            'summary_file': summary_path,
            'parameters': {
                'min_features': min_features,
                'max_features': max_features,
                'max_mt_percent': max_mt_percent,
                'n_variable_features': n_variable_features,
                'n_pcs': n_pcs,
                'clustering_resolution': clustering_resolution,
                'regress_out_vars': regress_out_vars,
                'use_sct_transform': use_sct_transform
            }
        }

        
def main():
    """
    Main function to run comprehensive QC analysis
    """
    # Initialize QC analyzer
    qc_analyzer = ComprehensiveSingleCellQC()
    
    # Set input and output paths
    input_h5ad = '/Users/xuz3/single_cell_python/Cux2_Atf4_paper_code/Atf4/data/processed/13_desc_v03.h5ad'
    output_dir = '/Users/xuz3/single_cell_python/Cux2_Atf4_paper_code/Atf4/qc_analysis_results'
    
    # Run comprehensive QC analysis
    try:
        results = qc_analyzer.generate_comprehensive_qc_report_v2(input_h5ad, output_dir)
        
        print("\\n🎉 QC Analysis completed successfully!")
        print("\\n📋 Generated files:")

        
    except Exception as e:
        print(f"❌ Error during QC analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()