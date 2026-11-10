# 3. Creating the Seurat object

The PBMC tutorial does this with:

```r
pbmc.data <- Read10X(data.dir = ".../filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

That doesn't work directly here, for two reasons:

1. **File format.** `Read10X()` expects files named exactly `matrix.mtx.gz`, `barcodes.tsv.gz`
   and `features.tsv.gz`. The GEO files have other names and the gene file has a header row.
   `ReadMtx()` is the equivalent that takes each file path separately.
2. **Size.** The PBMC dataset has 2,700 cells; this one has 312,928 cells and ~693 million
   non-zero values. Reading it into memory with `ReadMtx()` needs more than 16 GB of RAM.

## Solution: BPCells

BPCells converts the matrix once into an on-disk format that Seurat reads without loading
everything into memory. The full code is in [`scripts/01_load_data.R`](../scripts/01_load_data.R).

```r
d <- "data_raw/GSE136831"
genes <- read.delim(file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz"))[[2]]
cells <- readLines(file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz"))

# One-time conversion (10-30 min), writes data_raw/GSE136831/bpcells_counts/
counts <- import_matrix_market(file.path(d, "GSE136831_RawCounts_Sparse.mtx.gz"),
                               outdir = file.path(d, "bpcells_counts"),
                               row_names = make.unique(genes),
                               col_names = cells)

ipf <- CreateSeuratObject(counts = counts, project = "IPF")
```

- `make.unique()` is needed because some gene symbols occur more than once.
- `min.cells` / `min.features` are left out because they would read the full matrix into RAM;
  low-quality cells are removed in the QC step instead.

For comparison, the in-memory version (only for machines with plenty of RAM):

```r
ipf.data <- ReadMtx(mtx      = file.path(d, "GSE136831_RawCounts_Sparse.mtx.gz"),
                    cells    = file.path(d, "GSE136831_AllCells.cellBarcodes.txt.gz"),
                    features = file.path(d, "GSE136831_AllCells.GeneIDs.txt.gz"),
                    feature.column = 2, skip.feature = 1)
```

## Getting the object back in a new session

**Rebuild from the converted counts** (seconds):

```r
counts <- open_matrix_dir("data_raw/GSE136831/bpcells_counts")
ipf <- CreateSeuratObject(counts = counts, project = "IPF")
```

**Save and reload** (keeps QC, normalization, clusters, etc.):

```r
saveRDS(ipf, "ipf.rds")
ipf <- readRDS("ipf.rds")
```

With BPCells the `.rds` does not contain the counts, only a pointer to
`data_raw/GSE136831/bpcells_counts/`. Don't move or delete that folder. Once it exists, the original
`.mtx.gz` is no longer needed.
