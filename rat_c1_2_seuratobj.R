library(Seurat)
library(Matrix)

##Reading in data for Rat_C1_2
data_dir <- "X:/Chudasama49/Cordelia/Single Cell/Rat Unfiltered Matrixes/ec26219a-bf51-4ad0-af4e-042f545ebbe1_unfiltered_matrices/output_combined/Rat_C1_2/DGE_unfiltered"

setwd(data_dir)

expression.matrix <- ReadMtx(
  mtx = "count_matrix.mtx.gz",
  features = "all_genes.csv.gz",
  feature.column =1,
  cells = "cell_metadata.csv.gz",
  skip.cell = 1,
  skip.feature = 1,
  mtx.transpose = TRUE
)

seurat_object <-CreateSeuratObject(counts = expression.matrix)
