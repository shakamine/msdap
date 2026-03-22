make_model_term_dataset = function() {
  samples = tibble::tibble(
    sample_id = c("S1", "S2", "S3", "S4"),
    shortname = c("S1", "S2", "S3", "S4"),
    exclude = c(FALSE, FALSE, FALSE, FALSE),
    group = c("A", "A", "B", "B"),
    batch = c("b1", "b2", "b1", "b3"),
    donor = c("d1", "d2", "d3", "d4"),
    pair_id = c("p1", "p2", "p1", "p2"),
    score = c(0.1, 0.2, 0.3, 0.4)
  )

  proteins = tibble::tibble(protein_id = c("prot1", "prot2", "prot3"))
  peptides = tibble::tribble(
    ~sample_id, ~protein_id, ~peptide_id, ~sequence_plain, ~sequence_modified, ~detect, ~intensity,
    "S1", "prot1", "pep1", "PEP1", "PEP1", TRUE, 10,
    "S2", "prot1", "pep1", "PEP1", "PEP1", TRUE, 11,
    "S3", "prot1", "pep1", "PEP1", "PEP1", TRUE, 12,
    "S4", "prot1", "pep1", "PEP1", "PEP1", TRUE, 13,
    "S1", "prot2", "pep2", "PEP2", "PEP2", TRUE, 9,
    "S2", "prot2", "pep2", "PEP2", "PEP2", TRUE, 10,
    "S3", "prot2", "pep2", "PEP2", "PEP2", TRUE, 11,
    "S4", "prot2", "pep2", "PEP2", "PEP2", TRUE, 12,
    "S1", "prot3", "pep3", "PEP3", "PEP3", TRUE, 8,
    "S2", "prot3", "pep3", "PEP3", "PEP3", TRUE, 9,
    "S3", "prot3", "pep3", "PEP3", "PEP3", TRUE, 10,
    "S4", "prot3", "pep3", "PEP3", "PEP3", TRUE, 11
  )

  list(
    name = "model_terms_test",
    samples = samples,
    proteins = proteins,
    peptides = peptides,
    contrasts = list()
  )
}


make_limma_test_eset = function() {
  set.seed(123)
  mat = matrix(rnorm(12 * 4, mean = 10, sd = 1), nrow = 12, ncol = 4)
  mat[1:6, 3:4] = mat[1:6, 3:4] + 1
  rownames(mat) = paste0("prot", seq_len(nrow(mat)))
  colnames(mat) = c("S1", "S2", "S3", "S4")

  eset = Biobase::ExpressionSet(mat)
  Biobase::pData(eset) = data.frame(
    sample_id = colnames(mat),
    condition = c(1, 1, 2, 2),
    batch = c("b1", "b2", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  Biobase::fData(eset) = data.frame(
    protein_id = rownames(mat),
    npep = as.integer(rep(1:4, length.out = nrow(mat))),
    stringsAsFactors = FALSE,
    row.names = rownames(mat)
  )

  list(
    eset = eset,
    model_matrix = stats::model.matrix(~ condition + batch, data = Biobase::pData(eset)),
    block = c("p1", "p2", "p1", "p2")
  )
}


testthat::test_that("setup_contrasts stores fixed random and block model terms", {
  dataset = msdap::setup_contrasts(
    make_model_term_dataset(),
    contrast_list = list(c("A", "B")),
    fixed_variables = "batch",
    random_variables = "donor",
    block_variable = "pair_id"
  )

  contr = dataset$contrasts[[1]]
  testthat::expect_equal(contr$fixed_variables, "batch")
  testthat::expect_equal(contr$random_variables, "donor")
  testthat::expect_equal(contr$block_variable, "pair_id")
  testthat::expect_equal(contr$block_vector, c("p1", "p2", "p1", "p2"))
  testthat::expect_equal(contr$colname_additional_variables, "batch")
  testthat::expect_true("condition" %in% colnames(contr$model_matrix))
  testthat::expect_true(any(grepl("^batch", colnames(contr$model_matrix))))
  testthat::expect_false(any(grepl("donor", colnames(contr$model_matrix))))
  testthat::expect_true("donor" %in% colnames(contr$sample_table))
})


testthat::test_that("setup_contrasts rejects continuous numeric block variables", {
  testthat::expect_error(
    msdap::setup_contrasts(
      make_model_term_dataset(),
      contrast_list = list(c("A", "B")),
      block_variable = "score"
    )
  )
})


testthat::test_that("de_ebayes and de_deqms accept limma block vectors", {
  test_input = make_limma_test_eset()

  res_ebayes = msdap::de_ebayes(
    eset = test_input$eset,
    model_matrix = test_input$model_matrix,
    model_matrix_result_prop = "condition",
    limma_block_vector = test_input$block
  )
  testthat::expect_s3_class(res_ebayes, "tbl_df")
  testthat::expect_equal(nrow(res_ebayes), 12)

  res_deqms = msdap::de_deqms(
    eset = test_input$eset,
    model_matrix = test_input$model_matrix,
    model_matrix_result_prop = "condition",
    limma_block_vector = test_input$block
  )
  testthat::expect_s3_class(res_deqms, "tbl_df")
  testthat::expect_equal(nrow(res_deqms), 12)
})


testthat::test_that("dea forwards only random variables to plugin algorithms", {
  dataset = msdap::setup_contrasts(
    make_model_term_dataset(),
    contrast_list = list(c("A", "B")),
    fixed_variables = "batch",
    random_variables = "donor"
  )

  capture_env = new.env(parent = emptyenv())
  spy_terms = function(peptides, samples, eset_peptides, eset_proteins, model_matrix, model_matrix_result_prop, random_variables, dataset_name) {
    capture_env$random_variables = random_variables
    capture_env$model_columns = colnames(model_matrix)
    tibble::tibble(
      protein_id = unique(peptides$protein_id),
      pvalue = 1,
      qvalue = 1,
      foldchange.log2 = 0,
      dea_algorithm = "spy_terms"
    )
  }

  assign("spy_terms", spy_terms, envir = .GlobalEnv)
  on.exit(rm("spy_terms", envir = .GlobalEnv), add = TRUE)

  result = msdap::dea(dataset, dea_algorithm = "spy_terms", rollup_algorithm = "sum")
  testthat::expect_true("de_proteins" %in% names(result))
  testthat::expect_equal(capture_env$random_variables, "donor")
  testthat::expect_false(any(grepl("donor", capture_env$model_columns)))
  testthat::expect_true(any(grepl("^batch", capture_env$model_columns)))
})
