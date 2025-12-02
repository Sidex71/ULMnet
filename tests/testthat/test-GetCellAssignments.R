test_that("GetCellAssignments successfully classifies cells based on scores", {
  ###### Load test data from ULM package
  data("int_singData", package = "ULMnet") ## load scRNAseq data
  data("int_signature", package = "ULMnet") ## load signatures
  
  ##Compute signature scores on 200 cells
  my_scores <- GetCellScores(
    seurat_obj = int_singData[, 1:200],
    signatures = int_signature,
    assay = "RNA",
    layer = "data"
  )
  
  #######Run function for cell assignments
  my_ass <- GetCellAssignments(score_data = my_scores, cut_off = 1, p_val = 0.05)
  
  ## Basic checks on expected outputs
  expect_s3_class(my_ass, "data.frame")               ## output should be a data frame
  expect_true(all(c("barcode", "count_ulm", "celltype_ulm",
                    "avg_pvalue", "avg_score") %in% colnames(my_ass)))  # required columns
  expect_true(all(my_ass$avg_pvalue <= 0.05))         ## default p_val filter should apply during classification
  expect_true(all(my_ass$avg_score > 1))              ## default cut_off filter should apply during classification
  expect_true(all(!duplicated(my_ass$barcode)))       ## each barcode should appear once after classification
})
