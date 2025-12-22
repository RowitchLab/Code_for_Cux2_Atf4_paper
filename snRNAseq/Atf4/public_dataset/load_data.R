library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(loomR)

# Load Loom file
# Replace 'your_file.loom' with your actual file path
loom_file <- "~/single_cell_python/data/linnarsson_adolescent.loom"
# seurat_obj <- LoadH5Seurat(loom_file)
# Alternative method if above doesn't work:
# Connect to loom file and convert to Seurat
lfile <- Connect(filename = loom_file, mode = "r")
seurat_obj <- as.Seurat(lfile)

# View basic information about the object
print(seurat_obj)
print(dim(seurat_obj))