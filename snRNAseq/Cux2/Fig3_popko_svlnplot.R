
# Set working directory
setwd("/Users/xuz3/Library/CloudStorage/OneDrive-Cedars-SinaiHealthSystem/Documents/Code/work_related/Bioinfo/single_cell_python/Cux2_Atf4_paper_code/Cux2")
source('../utils/splitviolin_plot_base.R')
# Load Seurat object
seurat_object <- readRDS('./data/Popko_EN-L2-3.rds')

# Subset Seurat object
seurat_subset <- subset(seurat_object, subset= age %in% c("W5/7", "W17","W27/29","W41/44"))


# Add module scores for DDR, UPR, ISR, and AAT
seurat_subset <- add_module_score(seurat_subset, '../utils/DDR_list_final.csv', 'DDR')
seurat_subset <- add_module_score(seurat_subset, '../utils/IFN_sorted.csv', 'IFN')
seurat_subset <- add_module_score(seurat_subset, '../utils/UPR_list.csv', 'UPR')
seurat_subset <- add_module_score(seurat_subset, '../utils/ISR_list.csv', 'ISR')
# seurat_subset <- add_module_score(seurat_subset, '../utils/AAT.csv', 'AAT')
# seurat_subset <- add_module_score(seurat_subset, '../utils/PACT.csv', 'PACT')
# seurat_subset <- add_module_score(seurat_subset, '../utils/NRF2.csv', 'NRF2')

# Prepare metadata for plotting
metadata <- seurat_subset@meta.data
plot_data <- metadata[c( 'DDR1','IFN1', 'age', 'Condition','UPR1','ISR1')]

# List of genes to plot
gl_to_plot <- c('DDR1','UPR1','ISR1','IFN1')#,

# Generate plots
# for (tglname in gl_to_plot) {
#   plot_split_violin(plot_data, tglname,'age','Condition','_Popko_ENL2-3',fsizeh=4,fsizew=8)

# }
genes_to_plot <- c('Atf4','Cux2','Hdac9','Actb','Gpx4','2900097C17Rik','Apex1','Rad23b','Trp53bp1', 'Xrcc6')
#=================
# genes_to_plot <-c('2900097C17Rik','Cux2', 'Atf4', 'Actb', 'Hspa1a', 'Gpx4', 
#                 'Apex1', 'Parp1', 'Rad23b', 'Hdac9', 'Atm', 'Trp53bp1', 'Xrcc6','Xrcc5')
# genes_to_plot <- c('Rpa3','2900097C17Rik', 'Cux2', 'Atf4')
plot_data2 <- FetchData(seurat_subset, vars = c(genes_to_plot))
plot_data3 <- cbind(plot_data,plot_data2)

for (tglname in genes_to_plot) {
  plot_split_violin(plot_data3, tglname,'age','Condition','_Popko_ENL2-3',fsizeh=4,fsizew=8)
}


# genes_to_plot <- c('DDR1','IFN1','Rpa3','2900097C17Rik', 'Cux2', 'Atf4','Xrcc5')
for (tglname in genes_to_plot) {
  summarize_group_stats(plot_data3, tglname,'age','Condition','Popko_ENL2-3')
}
