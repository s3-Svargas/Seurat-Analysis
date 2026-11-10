# 2. Downloading the data from GEO

Dataset: **GSE136831**, the IPF Cell Atlas (IPF, COPD and control lungs).

## Metadata vs. counts

`getGEO()` only downloads the sample metadata. Single-cell count matrices are stored as
**supplementary files**, which are downloaded with `getGEOSuppFiles()`.

```r
gse_id <- "GSE136831"

# List the supplementary files without downloading
getGEOSuppFiles(gse_id, fetch_files = FALSE)
```

| File | Size | Contents | Needed? |
|---|---|---|---|
| `GSE136831_RawCounts_Sparse.mtx.gz` | 2.15 GB | raw UMI counts, all cells (Matrix Market) | **yes** |
| `GSE136831_AllCells.cellBarcodes.txt.gz` | 1.4 MB | 312,928 cell barcodes, no header | **yes** |
| `GSE136831_AllCells.GeneIDs.txt.gz` | 0.3 MB | 45,947 genes: Ensembl ID + gene symbol, with header | **yes** |
| `GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz` | 4.3 MB | sample, disease and cell type per cell | yes (annotation) |
| `GSE136831_RAW.tar` | 17.7 GB | per-sample `.dgecounts.rds.gz` files | no, same counts as the `.mtx` |

## Download

```r
options(timeout = 3600)   # the default 60 s timeout cuts off large downloads
getGEOSuppFiles(gse_id, baseDir = "data_raw",
                filter_regex = "RawCounts_Sparse|AllCells")
```

The `.mtx.gz` can also be downloaded in a browser from the GEO page and moved into
`data_raw/GSE136831/`.

Final layout:

```
data_raw/GSE136831/
├── GSE136831_AllCells.cellBarcodes.txt.gz
├── GSE136831_AllCells.GeneIDs.txt.gz
├── GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz
└── GSE136831_RawCounts_Sparse.mtx.gz
```

## Recognizing count matrices on GEO

File names that usually mean raw counts: `matrix.mtx` + `barcodes.tsv` + `features.tsv`/`genes.tsv`,
`.h5`, or names containing `counts`, `UMI`, `dge`, `raw`, `filtered`. Names with `normalized`,
`log`, `TPM`, `FPKM` or `CPM` are already transformed and are not suitable as Seurat input.
