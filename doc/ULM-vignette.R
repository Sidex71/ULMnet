## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ULM)

## -----------------------------------------------------------------------------
data("int_singData")
int_singData

## -----------------------------------------------------------------------------
head(int_singData@meta.data)

## ----fig.width=7, fig.height=6, fig.align='center'----------------------------
library(Seurat)
DimPlot(int_singData, reduction="umap", group.by="Cell_Type", label=TRUE)

## -----------------------------------------------------------------------------
set.seed(101324)
int_sig <- GetSignature(int_singData, ident_col = int_singData$Cell_Type, n = 100)

## -----------------------------------------------------------------------------
head(int_sig, 10)

## -----------------------------------------------------------------------------
my_scores <- GetCellScores(seurat_obj = int_singData, signatures = int_sig, assay = 'RNA', slot = 'data')

## -----------------------------------------------------------------------------
head(my_scores, 10)

## -----------------------------------------------------------------------------
my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1, p_val = 0.05)

## -----------------------------------------------------------------------------
head(my_ass, 10)

## -----------------------------------------------------------------------------
head(int_singData@meta.data)

## -----------------------------------------------------------------------------
int_singData <- AddMetaObject(int_singData, cell_class_df = my_ass)

## -----------------------------------------------------------------------------
head(int_singData@meta.data)

## -----------------------------------------------------------------------------
my_mult <- GetMultiplet(int_singData, minCells = 2)

## -----------------------------------------------------------------------------
multSummary <- my_mult$multSummary
multSummary

## -----------------------------------------------------------------------------
sum(multSummary$frequency)

## -----------------------------------------------------------------------------
multObj <- my_mult$multObj
multObj

## -----------------------------------------------------------------------------
table(multObj$count_ulm)

## -----------------------------------------------------------------------------
my_mult_filt <- FilterMultiplet(int_singData, minCells = 2, minFreq = 10)

## -----------------------------------------------------------------------------
multSummaryFilt <- my_mult_filt$multSummaryFilt
multSummaryFilt

## -----------------------------------------------------------------------------
multObjFilt <- my_mult_filt$multObjFilt
multObjFilt

## -----------------------------------------------------------------------------
my_node_df <- GetNodeDF(mat = multSummaryFilt)

## -----------------------------------------------------------------------------
my_node_df

## ----fig.width=8, fig.height=6, fig.align='left'------------------------------
PlotNetwork(my_node_df)

## -----------------------------------------------------------------------------
data("int_singData")
ref_obj <- int_singData
ref_obj

## -----------------------------------------------------------------------------
data("int_multData")
query_obj <- int_multData
query_obj

## -----------------------------------------------------------------------------
set.seed(101324)
ref_sig <- GetSignature(ref_obj, ident_col = ref_obj$Cell_Type)
head(ref_sig, 10)

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
my_scores <- GetCellScores(seurat_obj = query_obj, signatures = ref_sig, assay = 'RNA', slot = 'data')

## -----------------------------------------------------------------------------
my_ass <- GetCellAssignments(score_data = my_scores)

## -----------------------------------------------------------------------------
lab_query_obj <- AddMetaObject(query_obj, cell_class_df = my_ass)

## -----------------------------------------------------------------------------
query_mult <- GetMultiplet(lab_query_obj)

## -----------------------------------------------------------------------------
query_multSummary <- query_mult$multSummary
query_multSummary

## -----------------------------------------------------------------------------
query_multObj <- query_mult$multObj
query_multObj
table(query_multObj$count_ulm)

## -----------------------------------------------------------------------------
query_mult_filt <- FilterMultiplet(lab_query_obj, minFreq = 7)

## -----------------------------------------------------------------------------
query_multSummaryFilt <- query_mult_filt$multSummaryFilt
query_multSummaryFilt

## -----------------------------------------------------------------------------
query_network_df <- GetNodeDF(mat = query_multSummaryFilt)

## ----fig.width=15, fig.height=12, fig.align='left'----------------------------
PlotNetwork(query_network_df)

## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


