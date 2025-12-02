test_that("GetSignature generates signatures correctly", {
  # load an example scRNAseq data provided with the package
  data("int_singData", package = "ULMnet")
  
  # run the function with a specified annotation column
  expect_message(
    int_sig <- GetSignature(
      seurat_obj = int_singData[, 1:500], ###using 500 cells
      ident_col = int_singData$Cell_Type,
      n = 100,
      p_val = 0.05
    ),
    regexp = "using the specified seurat ident"
  )
  
  # structure tests
  expect_s3_class(int_sig, "data.frame")
  expect_true(all(c("source", "target", "mor") %in% colnames(int_sig)))
  
  # check the content of each column
  expect_true(all(int_sig$mor == 1))
  expect_true(all(!is.na(int_sig$source)))
  expect_true(all(!is.na(int_sig$target)))
  
  # check that the number of genes in each cell type signature is ≤ n
  cluster_counts <- table(int_sig$source)
  expect_true(all(cluster_counts <= 100))
  
})

