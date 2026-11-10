# Seurat-Analysis
Utilizing Seurat for data analysis of scRNA-seq

A step-by-step tutorial for taking **any single-cell RNA-seq dataset from GEO** and turning it into a
Seurat object, ready for the [Seurat PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
workflow (QC → normalization → clustering → markers). It covers every common GEO file type.

The worked example is the IPF Cell Atlas, [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831)
(312,928 cells from IPF, COPD and control lungs). Details for that dataset are in [`docs/`](docs/).

---

## Contents

1. [Set up R](#1-set-up-r)
2. [Find the count files on GEO](#2-find-the-count-files-on-geo)
3. [Download](#3-download)
4. [Read the counts — by file type](#4-read-the-counts--by-file-type)
5. [Large datasets: BPCells](#5-large-datasets-bpcells)
6. [Create the Seurat object](#6-create-the-seurat-object)
7. [Add the GEO metadata](#7-add-the-geo-metadata)
8. [Save and reload](#8-save-and-reload)
9. [Next: the PBMC tutorial steps](#9-next-the-pbmc-tutorial-steps)

---

## 1. Set up R

```r
install.packages(c("Seurat", "patchwork", "dplyr", "ggplot2", "BiocManager", "data.table"))
BiocManager::install("GEOquery")
install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", getOption("repos")))
```

Optional, only if your dataset uses these formats:

```r
BiocManager::install("zellkonverter")    # .h5ad (AnnData / scanpy)
BiocManager::install("LoomExperiment")   # .loom
```

This repo's [`.Rprofile`](.Rprofile) loads the main packages automatically whenever R starts in this
folder. See [docs/01_setup.md](docs/01_setup.md) for how that works and how to update packages safely.

> Write code in an R script (File → New File → R Script) and run it with **Cmd+Enter**. Code pasted into
> the console runs immediately, line by line.

## 2. Find the count files on GEO

`getGEO()` only downloads **sample metadata**. The count matrix is in the **supplementary files**.
List them without downloading:

```r
library(GEOquery)
gse_id <- "GSE136831"                          # ← your accession
files <- getGEOSuppFiles(gse_id, fetch_files = FALSE)
files$fname
```

Pick what to download by file name:

| File name looks like | What it is | Use it? |
|---|---|---|
| `matrix.mtx(.gz)` + `barcodes.tsv` + `features.tsv` / `genes.tsv` | 10x / Matrix Market counts | ✅ best choice |
| `*.mtx.gz` with other names (`RawCounts_Sparse.mtx.gz`, `cellBarcodes.txt`, `GeneIDs.txt`) | same, custom names | ✅ |
| `*.h5` (`filtered_feature_bc_matrix.h5`) | 10x HDF5 counts | ✅ |
| `*counts*`, `*UMI*`, `*dge*` as `.csv` / `.tsv` / `.txt(.gz)` | dense table of counts | ✅ |
| `*.rds` | R object (often a finished Seurat object) | ✅ but often huge and already processed |
| `*.h5ad` | AnnData (Python / scanpy) | ✅ |
| `*.loom` | loom file | ✅ |
| `*_RAW.tar` | archive of per-sample files | only if nothing above exists |
| `*normalized*`, `*log*`, `*TPM*`, `*FPKM*`, `*CPM*` | already transformed values | ❌ Seurat needs raw counts |
| `*metadata*`, `*annotation*`, `*celltype*` | per-cell info | ✅ alongside the counts (step 7) |

Prefer **filtered** over **raw** for 10x files: "raw" includes empty droplets.

## 3. Download

```r
options(timeout = 3600)   # the default 60 s timeout cuts off large files

# Download only the files you need (regex on file names)
getGEOSuppFiles(gse_id, baseDir = "data_raw",
                filter_regex = "matrix|barcodes|features|genes|metadata")
list.files(file.path("data_raw", gse_id))
```

- Files land in `data_raw/<GSE>/`, relative to `getwd()`.
- If the R download fails, download from the GEO web page in a browser and move the file into `data_raw/<GSE>/`.
- Check you have enough free disk space first; `_RAW.tar` files can be tens of GB.

For a `_RAW.tar`, unpack it:

```r
untar(file.path("data_raw", gse_id, paste0(gse_id, "_RAW.tar")),
      exdir = file.path("data_raw", gse_id, "raw"))
```

## 4. Read the counts — by file type

Every option below ends with `counts`: a **genes × cells** matrix with gene names as row names and
cell barcodes as column names. Set the folder first:

```r
d <- file.path("data_raw", gse_id)
```

### A. 10x folder with standard names

The folder contains exactly `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` (or `genes.tsv.gz`):

```r
counts <- Read10X(data.dir = d)
```

### B. Matrix Market with custom names or prefixes

GEO usually renames files (`GSM123_matrix.mtx.gz`, `RawCounts_Sparse.mtx.gz`, …), so `Read10X` can't
find them. Use `ReadMtx` and point at each file:

```r
counts <- ReadMtx(mtx      = file.path(d, "GSE136831_RawCounts_Sparse.mtx.gz"),
                  cells    = file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz"),
                  features = file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz"),
                  feature.column = 2,   # column holding gene symbols (1 = Ensembl ID)
                  skip.feature = 1,     # 1 if the gene file has a header row
                  skip.cell = 0)        # 1 if the barcode file has a header row
```

Peek at the files to choose these arguments:

```r
readLines(gzfile(file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz")), n = 3)
readLines(gzfile(file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz")), n = 3)
```

### C. 10x HDF5 (`.h5`)

```r
counts <- Read10X_h5(file.path(d, "filtered_feature_bc_matrix.h5"))
# Multimodal files (e.g. CITE-seq) return a list:
# counts <- counts$`Gene Expression`
```

### D. Dense table (`.csv`, `.tsv`, `.txt`)

```r
tab <- data.table::fread(file.path(d, "counts.csv.gz"))
counts <- as.matrix(tab[, -1])
rownames(counts) <- tab[[1]]
counts <- as(counts, "dgCMatrix")          # sparse: much less memory

counts[1:5, 1:5]                           # genes should be rows
# counts <- t(counts)                      # flip if cells are rows
```

### E. R object (`.rds`)

```r
obj <- readRDS(file.path(d, "data.rds"))
class(obj)
```

- `"Seurat"` → you already have an object: `ipf <- UpdateSeuratObject(obj)` and skip to step 7.
  It may already be normalized and clustered.
- `"dgCMatrix"` or `"matrix"` → `counts <- obj`.

### F. AnnData (`.h5ad`)

```r
sce <- zellkonverter::readH5AD(file.path(d, "data.h5ad"))
SummarizedExperiment::assayNames(sce)     # see which assay holds raw counts
ipf <- as.Seurat(sce, counts = "X", data = NULL)   # use the raw-count assay name
```

In many `.h5ad` files `X` is normalized and raw counts are in a layer named `counts`; check before
converting (step 4 check below).

### G. loom

```r
lx <- LoomExperiment::import(file.path(d, "data.loom"), type = "SingleCellLoomExperiment")
counts <- as(SummarizedExperiment::assay(lx, "matrix"), "dgCMatrix")
rownames(counts) <- SummarizedExperiment::rowData(lx)$Gene   # column names vary by file
colnames(counts) <- SummarizedExperiment::colData(lx)$CellID
```

### H. One set of files per sample (`_RAW.tar`)

Read each sample with the matching option above, make one object per sample, then merge:

```r
mtx_files <- list.files(file.path(d, "raw"), pattern = "matrix.mtx.gz$", full.names = TRUE)

objs <- lapply(mtx_files, function(m) {
  prefix <- sub("_?matrix.mtx.gz$", "", m)
  cnt <- ReadMtx(mtx = m,
                 cells    = paste0(prefix, "_barcodes.tsv.gz"),
                 features = paste0(prefix, "_features.tsv.gz"))
  CreateSeuratObject(cnt, project = basename(prefix), min.cells = 3, min.features = 200)
})

ipf <- merge(objs[[1]], objs[-1], add.cell.ids = sapply(objs, Project))
```

`add.cell.ids` keeps barcodes unique across samples. Adjust the file-name patterns to your files.

### Check: are these raw counts?

```r
x <- counts[, 1:100]
range(x)
all(x == round(x))   # TRUE  → raw counts, good
                     # FALSE → normalized; look for another file
```

## 5. Large datasets: BPCells

The PBMC tutorial has 2,700 cells. Datasets with more than ~100k cells may not fit in RAM: reading
GSE136831's 313k cells with `ReadMtx` needs more than 16 GB. **BPCells** converts the counts once into a fast
on-disk format, and Seurat reads from disk.

```r
library(BPCells)
bp_dir <- file.path(d, "bpcells_counts")
```

Pick the import that matches your file type (one-time; can take 10–30 min):

```r
# Matrix Market (.mtx / .mtx.gz)
genes <- read.delim(file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz"))[[2]]
cells <- readLines(file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz"))
counts <- import_matrix_market(file.path(d, "GSE136831_RawCounts_Sparse.mtx.gz"),
                               outdir = bp_dir,
                               row_names = make.unique(genes),
                               col_names = cells)

# 10x HDF5 (.h5)
counts <- open_matrix_10x_hdf5(file.path(d, "filtered_feature_bc_matrix.h5"))
counts <- write_matrix_dir(counts, dir = bp_dir)

# AnnData (.h5ad), raw counts in X
counts <- open_matrix_anndata_hdf5(file.path(d, "data.h5ad"))
counts <- write_matrix_dir(counts, dir = bp_dir)
```

In later sessions, open it in seconds instead of converting again:

```r
counts <- open_matrix_dir(bp_dir)
```

- `make.unique()` avoids errors from duplicated gene symbols.
- Keep `bpcells_counts/` where it is: Seurat objects built on it only store its location.
- Once it exists, the original download can be deleted to save space.

## 6. Create the Seurat object

Same as the PBMC tutorial's `CreateSeuratObject` step:

```r
ipf <- CreateSeuratObject(counts = counts, project = "IPF",
                          min.cells = 3, min.features = 200)
ipf
```

With BPCells counts, leave out `min.cells` and `min.features` (they load the whole matrix into
memory) and filter cells in the QC step instead:

```r
ipf <- CreateSeuratObject(counts = counts, project = "IPF")
```

## 7. Add the GEO metadata

Per-cell metadata (sample, disease, cell type) usually comes as a separate supplementary file:

```r
meta <- data.table::fread(file.path(d, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"),
                          data.table = FALSE)
head(meta)

rownames(meta) <- meta[[1]]                 # the column holding cell barcodes
meta <- meta[colnames(ipf), , drop = FALSE] # same order as the object
ipf <- AddMetaData(ipf, meta)
head(ipf[[]])
```

Per-sample information (patient, diagnosis, age) is in the GEO series matrix:

```r
gse <- getGEO(gse_id, GSEMatrix = TRUE)
pheno <- Biobase::pData(gse[[1]])
```

Match it to cells by sample ID before adding it with `AddMetaData`.

## 8. Save and reload

```r
saveRDS(ipf, "ipf.rds")      # end of a session or after each big step
ipf <- readRDS("ipf.rds")    # next session
```

Run R from the same folder each time so the relative paths (`data_raw/…`, `ipf.rds`) still work.

## 9. Next: the PBMC tutorial steps

From here the [PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) applies
directly, with `pbmc` replaced by your object:

```r
ipf[["percent.mt"]] <- PercentageFeatureSet(ipf, pattern = "^MT-")   # mouse: "^mt-"
VlnPlot(ipf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

then QC filtering → `NormalizeData` → `FindVariableFeatures` → `ScaleData` → `RunPCA` →
`FindNeighbors` / `FindClusters` → `RunUMAP` → `FindAllMarkers`.

For very large objects, see Seurat's [sketch-based analysis](https://satijalab.org/seurat/articles/seurat5_sketch_analysis)
for running these steps on a representative subset of cells.

---

## Progress on GSE136831

- [x] R environment set up — [docs/01_setup.md](docs/01_setup.md)
- [x] Data downloaded from GEO — [docs/02_data_download.md](docs/02_data_download.md)
- [x] Counts converted with BPCells, Seurat object created — [docs/03_create_seurat_object.md](docs/03_create_seurat_object.md), [`scripts/01_load_data.R`](scripts/01_load_data.R)
- [ ] QC and cell filtering
- [ ] Normalization, variable features, scaling
- [ ] PCA, clustering, UMAP
- [ ] Marker genes and cell type annotation

## Layout

```
.Rprofile                 auto-loads the analysis packages in interactive R sessions
scripts/01_load_data.R    builds the Seurat object for GSE136831
docs/                     notes for each step of the GSE136831 analysis
data_raw/                 GEO downloads (not tracked by git, see .gitignore)
```

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
