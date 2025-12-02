test_that("GetMultiplet correctly extracts predicted multiplets", {
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
  
  ################ run function to extract multiplets (expect warning about removing missing cells)
  
  expect_warning(
    my_mult <- GetMultiplet(seurat_obj = new_obj, minCells = 2),
    regexp = "Removing .* cells missing data"
  )
  ### Basic checks on expected outputs
  expect_type(my_mult, "list")                           ## output is a list
  expect_true(all(c("multSummary", "multObj") %in% names(my_mult))) ## the list has the expected 2 components
  
  ## checks on multiplet summary output
  expect_s3_class(my_mult$multSummary, "data.frame")     ## summary is a data frame
  expect_true(all(c("multipletType", "frequency") %in% colnames(my_mult$multSummary))) ##must contain multipletType and frequency columns
  
  ## checks on multiplet Seurat object
  expect_s4_class(my_mult$multObj, "Seurat")             # subset should still be a Seurat object
  expect_true(all(my_mult$multObj$count_ulm >= 2))       # all barcodes satisfy minCells filter
})
