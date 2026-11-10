# Build a Seurat object from the GSE136831 (IPF Cell Atlas) GEO files.
# Run with the repo root as the working directory. See docs/03_create_seurat_object.md.

library(Seurat)
library(BPCells)

d <- "data_raw/GSE136831"
bp_dir <- file.path(d, "bpcells_counts")

if (dir.exists(bp_dir)) {
  # Already converted: open the on-disk matrix
  counts <- open_matrix_dir(bp_dir)
} else {
  # One-time conversion of the Matrix Market file (10-30 min)
  genes <- read.delim(file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz"))[[2]]
  cells <- readLines(file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz"))
  counts <- import_matrix_market(file.path(d, "GSE136831_RawCounts_Sparse.mtx.gz"),
                                 outdir = bp_dir,
                                 row_names = make.unique(genes),
                                 col_names = cells)
}

ipf <- CreateSeuratObject(counts = counts, project = "IPF")
ipf

saveRDS(ipf, "ipf.rds")
