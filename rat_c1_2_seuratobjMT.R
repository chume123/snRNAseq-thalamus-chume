library(Seurat)
library(Matrix)


input_dir <- normalizePath("X:/Chudasama49/Cordelia/Single Cell/RatUnfilteredMatrixes_MT/output_combined/C1_2/DGE_unfiltered")

obj <- Seurat::ReadParseBio(input_dir)

head(obj)

colnames(obj) |> head()

rownames(obj)

#need to add meta data etc,
cell_meta <- read.csv(paste0(input_dir, "/cell_metadata.csv"), row.names = 1)


seurat_obj <- CreateSeuratObject(obj, names.field = 0, meta.data = cell_meta)
head(seurat_obj)

#setting initial identity class to Rat_C1_2
seurat_obj@meta.data$orig.ident <- factor(rep("rat_c1_2", nrow(seurat_obj@meta.data)))
Idents(seurat_obj) <- seurat_obj@meta.data$orig.ident
