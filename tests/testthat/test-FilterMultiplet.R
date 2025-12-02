

test_that("FilterMultiplet filters multiplets correctly", {
  ###### Load test data from ULM package
  data("int_multData", package = "ULMnet")   ## load scRNAseq data
  data("int_signature", package = "ULMnet")  ## load signatures
  
  ######### prepare data: run upstream steps 
  ## Compute signature scores on 300 cells
  my_scores <- GetCellScores(
    seurat_obj = int_multData[, 1:300],
    signatures = int_signature,
    assay = "RNA",
    layer = "data"  
  )
  
  ## Assign cell types
  my_ass <- GetCellAssignments(score_data = my_scores)
  
  ## Add assignments to metadata
  new_obj <- AddMetaObject(seurat_obj = int_multData[, 1:300], cell_class_df = my_ass)

  #########################Run function to filter multiplets (expect a warning about removed cells)
  expect_warning(
    my_mult_filt <- FilterMultiplet(seurat_obj = new_obj, minCells = 2, minFreq = 2),
    regexp = "Removing .* cells missing data"
  )
  ####Basic checks on expected outputs
  expect_type(my_mult_filt, "list")  #                            #must be a list
  expect_named(my_mult_filt, c("multSummaryFilt", "multObjFilt")) ##list must have 2 components
  
  ## checks on filtered multiplet summary output
  expect_s3_class(my_mult_filt$multSummaryFilt, "data.frame")          ###must be a data frame
  expect_true(all(c("multipletType", "frequency") %in% colnames(my_mult_filt$multSummaryFilt))) ##must contain multipletType and frequency columns
  
  ## checks on multiplet Seurat object
  expect_s4_class(my_mult_filt$multObjFilt, "Seurat")             ##must be a Seurat object
  expect_true(all(my_mult_filt$multSummaryFilt$frequency >= 2))   ## all barcodes satisfy minCells filter
})
