# Auto-load packages in interactive sessions only (Rscript runs are unaffected).
# To update packages, start R with `R --no-init-file` so nothing is preloaded.
if (interactive()) {
  suppressMessages({
    library(Seurat)
    library(patchwork)
    library(dplyr)
    library(ggplot2)
    library(BiocManager)
    library(GEOquery)
    library(BPCells)
  })
}
