test_that("AddMetaObject correctly adds cell assignments to metadata", {
  ###### Load test data from ULM package
  data("int_singData", package = "ULMnet") ## load scRNAseq data
  data("int_signature", package = "ULMnet") ## load signatures
  
  ######### prepare data: run upstream steps 
  ## Compute signature scores on 200 cells
  my_scores <- GetCellScores(
    seurat_obj = int_singData[, 1:200],
    signatures = int_signature,
    assay = "RNA",
    layer = "data" 
  )
  
  ## Assign cell types
  my_ass <- GetCellAssignments(score_data = my_scores)
  
  ################ run function to add assignments to metadata
  new_obj <- AddMetaObject(seurat_obj = int_singData[, 1:200], cell_class_df = my_ass)
  
  ## Basic checks on expected outputs
  expect_s4_class(new_obj, "Seurat")   # object must still be a Seurat object
  
  expect_true(all(c("celltype_ulm", "count_ulm") %in% colnames(new_obj@meta.data))) ## new columns must be added
  
  expect_equal(
    sum(!is.na(new_obj$celltype_ulm)),
    nrow(my_ass)   ## there must be same number of assignments as input df
  )
  
  expect_true(all(my_ass$celltype_ulm %in% new_obj$celltype_ulm)) ## cell labels should match
})
