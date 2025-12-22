#%%
library(Seurat)
# setwd('/home/john/Code/Bioinfo/single_cell_R/Cux2-Atf4/Seurat_related/load_data')
setwd('/Users/xuz3/Code/work_related/Bioinfo/single_cell_python/Cux2_Atf4_paper_code/Cux2')
library(ggplot2)
# library(introdataviz)
seurat_object <- readRDS('./data/Lucas.rds')    # Close the PDF device



#%%

seurat_subset <- subset(seurat_object, subset= cell_type %in% c("EN-L2-3-A", "EN-L2-3-B"))

my_cols <- c("#7A7A7A", "#DC9F41") 

genes_to_plot <- c('HDAC9','USP11', 'ACTB','GPX4','HSPA1A','ATM', 'APEX1', 'XRCC6','PARP1', 'RAD23B',  'TP53BP1' )
for (tgene in genes_to_plot){
    violin_plot<- VlnPlot(
    seurat_subset,
    features = tgene,
    split.by = "diagnosis",
    pt.size=0.02,
    split.plot=TRUE
    )+ theme_minimal() +
    theme(
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    axis.line = element_line(color = "black")  # Adds black axis lines
  )+
    scale_fill_manual(values = my_cols) 

    ggsave(
    filename = file_name <- paste0("figures/TSplit_ViolinPlot_", tgene, ".pdf"),
    plot = violin_plot,
    width = 4,
    height = 4
    )
}
