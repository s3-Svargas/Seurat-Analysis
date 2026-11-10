# 1. R environment setup

## Bioconductor

GEOquery is a Bioconductor package, so it is installed through BiocManager.

```r
install.packages("BiocManager")
BiocManager::version()          # 3.22 for R 4.5
BiocManager::install("GEOquery")
```

To check whether a package is installed:

```r
requireNamespace("BiocManager", quietly = TRUE)   # TRUE if installed
packageVersion("BiocManager")
```

## BPCells

BPCells stores large count matrices on disk instead of in RAM. It is installed from r-universe:

```r
install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", getOption("repos")))
```

## Auto-loading packages with `.Rprofile`

R runs a file named `.Rprofile` from the working directory every time it starts. The
[`.Rprofile`](../.Rprofile) in this repo loads the analysis packages, so `library()` calls are
not needed each session:

```r
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
```

Notes:

- It only applies when R starts **in this folder**. R sessions elsewhere load nothing extra.
- `interactive()` means `Rscript` runs and knitting are unaffected, so scripts still need their own
  `library()` calls (the scripts in `scripts/` include them).
- Loaded packages appear in the **Packages** pane and in `search()`, not in the Environment pane,
  which only shows your own objects.
- To edit the list: `file.edit(".Rprofile")`, then restart R.

## Updating packages

A package that is already loaded can't be replaced, which gives errors like
`namespace 'rlang' 1.1.6 is already loaded, but >= 1.3.0 is required`. Because `.Rprofile`
loads packages at startup, update packages from a session that skips it:

```sh
R --no-init-file
```

then run `update.packages()` or `install.packages(...)`, and restart R normally.

## Running code

Paste code into an R script (File > New File > R Script) rather than the console. The console
runs every pasted line immediately; in a script you run lines with **Cmd+Enter** (current line
or selection) or **Cmd+Shift+Enter** (whole file).
