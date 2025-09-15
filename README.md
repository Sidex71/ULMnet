
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ULM

<!-- badges: start -->
<!-- badges: end -->

An R package for the prediction and identification of multiplets from
scRNAseq datasets to infer physical cell-cell interaction networks.
Multiplets occur naturally in conventional scRNAseq due to incomplete
dissociation during library preparation. They represent cells which are
physically connected and interacting in tissues which become sequenced
together as they remain unseperated. ULM utilizes a signature-based
approach where univariate linear models are fitted over each barcode in
a scRNAseq data to assign signature scores. Barcodes are then classified
as singlets or multiplets based on their signature scores. Multiplets
are those barcodes or cells that are enriched in two or more cell
type-specific gene signatures.

<img src="man/figures/algorithm.jpeg" style="width:100.0%" />

## Installation

You can install the development version of ULM from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

devtools::install_github("Sidex71/ULM", build_vignettes = T)
```

A full guide on ULM workflow can be found on the package website
<https://sidex71.github.io/ULM/articles/ULM-vignette.html> or on the
vignettes page by running

``` r
browseVignettes('ULM') 
#> No vignettes found by browseVignettes("ULM")
```

## Dependencies

This package depends on the following packages: Seurat (\>= 3.0.0),
decoupleR, tidyverse, dplyr, stringr, igraph, ggraph, tidygraph,
ggplot2, magrittr, tibble All the dependencies might be installed
alongside the package.

## Example

This is a quick example which shows how to infer physical cell-cell
interaction network from a scRNAseq data:

``` r
library(ULM)
####load dataset
data("int_singData")  ##int_singData is a preprocessed scRNAseq seurat object with a Cell_Type column containing cell annotations.
##generate signatures
set.seed(101324)
int_sig <- ULM::GetSignature(int_singData, ident_col = int_singData$Cell_Type, n = 100)
#> using the specified seurat ident to generate signatures
#> Calculating cluster Progenitor early
#> Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
#> ℹ Please use the `layer` argument instead.
#> ℹ The deprecated feature was likely used in the Seurat package.
#>   Please report the issue at <https://github.com/satijalab/seurat/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> Warning: `PackageCheck()` was deprecated in SeuratObject 5.0.0.
#> ℹ Please use `rlang::check_installed()` instead.
#> ℹ The deprecated feature was likely used in the Seurat package.
#>   Please report the issue at <https://github.com/satijalab/seurat/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
#> Calculating cluster Progenitor late-1
#> Calculating cluster Transit amplifying
#> Calculating cluster Progenitor late-2
#> Calculating cluster Goblet
#> Calculating cluster Stem
#> Calculating cluster Enterocyte
#> Calculating cluster Paneth
#> Calculating cluster Enteroendocrine
#> Calculating cluster Tuft

###score cells and assign labels
my_scores <- GetCellScores(seurat_obj = int_singData, signatures = int_sig, assay = 'RNA', layer  = 'data')
my_assign <- GetCellAssignments(score_data = my_scores, cut_off = 1, p_val = 0.05)
int_singData <- AddMetaObject(int_singData, cell_class_df = my_assign)

#####filter multiplets
my_mult_filt <- FilterMultiplet(int_singData, minCells = 2, minFreq = 10)
#> Warning: Removing 1929 cells missing data for vars requested
multSummaryFilt <- my_mult_filt$multSummaryFilt

###plot network
my_node_df <- GetNodeDF(mat = multSummaryFilt)
PlotNetwork(my_node_df)
```

<img src="man/figures/README-example-1.png" width="100%" />

A more comprehensive guide on the usage of the ULM package can be found
in the following link
<https://sidex71.github.io/ULM/articles/ULM-vignette.html>
