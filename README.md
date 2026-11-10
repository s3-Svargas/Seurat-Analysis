# Seurat-Analysis
Utilizing Seurat for data analysis of scRNA-seq

This repo follows the [Seurat PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
workflow, applied to the IPF Cell Atlas dataset from GEO
([GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831)): 312,928 cells from
idiopathic pulmonary fibrosis (IPF), COPD and control lungs.

## Progress

- [x] R environment set up (Bioconductor, auto-loaded packages) — [docs/01_setup.md](docs/01_setup.md)
- [x] Dataset downloaded from GEO — [docs/02_data_download.md](docs/02_data_download.md)
- [x] Count matrix converted with BPCells and loaded into a Seurat object — [docs/03_create_seurat_object.md](docs/03_create_seurat_object.md)
- [ ] QC and cell filtering
- [ ] Normalization, variable features, scaling
- [ ] PCA, clustering, UMAP
- [ ] Marker genes and cell type annotation

## Layout

```
.Rprofile                      auto-loads the analysis packages in interactive R sessions
scripts/01_load_data.R         builds the Seurat object from the GEO files
docs/                          step-by-step notes
data_raw/                      GEO downloads (not tracked by git, see .gitignore)
```

## Quick start

1. Open R with this folder as the working directory (e.g. open the folder as an RStudio project).
2. Put the GEO files in `data_raw/GSE136831/` (see [docs/02_data_download.md](docs/02_data_download.md)).
3. Run `scripts/01_load_data.R`.

## Software

| Package | Version |
|---|---|
| R | 4.5.1 |
| Seurat | 6.0.0 |
| BPCells | 0.3.1 |
| GEOquery | 2.78.0 |
| BiocManager | 1.30.27 (Bioconductor 3.22) |
| patchwork | 1.3.2 |
| dplyr | 1.2.1 |
| ggplot2 | 4.0.3 |
