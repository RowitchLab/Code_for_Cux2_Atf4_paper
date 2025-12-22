library(plot1cell)
library(Seurat)
library(ggplot2)
setwd("/Users/xuz3/single_cell_python/Cux2_Atf4_paper_code/Atf4/public_dataset")


seurat_object <- readRDS("./linnarsson_adolescent_filtered.rds")

print(seurat_object)
Idents(seurat_object) <- "TaxonomyRank4"
pdf("violin_Age_Atf4.pdf", width=5, height=4)
complex_vlnplot_single(seurat_object, feature = "Atf4", groups = "Age")
dev.off()


seurat_object <- readRDS("./linnarsson_dev_all_filtered.rds")

print(seurat_object)
Idents(seurat_object) <- "Age"
pdf("violin_Age_Atf4_dev.pdf", width=8, height=8)
complex_vlnplot_single(seurat_object, feature = "Atf4", groups = "Subclass",font.size=10)+ theme(axis.text.y = element_text(angle = 75, size=10, vjust = 1))
dev.off()


pdf("violin_Age_Atf4_dev2.pdf", width=8, height=8)
p <- complex_vlnplot_single(seurat_object, feature = "Atf4", groups = "Subclass",font.size=10)


p + ggplot2::theme(axis.text.y = ggplot2::element_text(vjust = 0.5, angle = 75))

print(p)
dev.off()
# # 3. Increase the space between the axis text and the plot
# p + ggplot2::theme(
#   axis.text.x = ggplot2::element_text(
#     vjust = 0.5,
#     margin = ggplot2::margin(t = 10) # Add a top margin of 10 points
#   )
# )

seurat_object@meta.data$Subclass <- factor(seurat_object@meta.data$Subclass,levels =c( "Cerebelllum glutamatergic",
                                          "Cerebelllum GABAergic",
                                        "Forebrain glutamatergic",
                                        "Forebrain GABAergic",
                                        "Midbrain glutamatergic",
                                        "Midbrain GABAergic",
                                        "Cortical or hippocampal glutamatergic"))
Idents(seurat_object) <- "Age"
pdf("violin_Age_Atf4_dev.pdf", width=8, height=8)
complex_vlnplot_single(seurat_object, feature = "Atf4", groups = "Subclass",font.size=10)+ theme(axis.text.y = element_text(angle = 75, size=10, vjust = 1))
dev.off()