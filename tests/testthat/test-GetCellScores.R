test_that("GetCellScores returns cell signature scores", {
  ### Load test data from the ULM package
  data("int_singData", package = "ULMnet") ## load scRNAseq data
  data("int_signature", package = "ULMnet") ## load signatures
  
  ### Run function on 200 cells using the normalised data layer
  scores <- GetCellScores(
    seurat_obj = int_singData[, 1:200],
    signatures = int_signature,
    assay = "RNA",
    layer = "data"
  )
  
  ### Structure checks
  expect_s3_class(scores, "data.frame")
  expect_true(all(c("barcode", "celltype", "score", "p_value", "statistic") %in% colnames(scores)))
  expect_gt(nrow(scores), 0)
  
  ### Sanity checks on values
  expect_true(all(is.finite(scores$score)))
  expect_true(all(is.finite(scores$p_value)))
  

  
})
