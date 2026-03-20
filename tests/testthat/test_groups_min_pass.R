testthat::context("groups_min_pass filtering")
msdap::enable_log(FALSE)

make_groups_min_pass_dataset = function() {
  sid = c("A1", "A2", "B1", "B2", "C1", "C2")
  samples = tibble::tibble(
    sample_id = sid,
    shortname = sid,
    group = c("A", "A", "B", "B", "C", "C"),
    exclude = FALSE
  )

  peptides = dplyr::bind_rows(
    tibble::tibble(
      protein_id = "p1",
      peptide_id = "pep_all",
      sequence_plain = "AAAAAA",
      sequence_modified = "AAAAAA",
      isdecoy = FALSE,
      detect = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
      rt = NA_real_,
      sample_id = sid,
      intensity = c(10, 10, 10, 10, 10, 10)
    ),
    tibble::tibble(
      protein_id = "p1",
      peptide_id = "pep_ab",
      sequence_plain = "BBBBBB",
      sequence_modified = "BBBBBB",
      isdecoy = FALSE,
      detect = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
      rt = NA_real_,
      sample_id = sid,
      intensity = c(20, 20, 20, 20, 20, 20)
    ),
    tibble::tibble(
      protein_id = "p1",
      peptide_id = "pep_a",
      sequence_plain = "CCCCCC",
      sequence_modified = "CCCCCC",
      isdecoy = FALSE,
      detect = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
      rt = NA_real_,
      sample_id = sid,
      intensity = c(30, 30, 30, 30, 30, 30)
    )
  )

  proteins = tibble::tibble(
    protein_id = "p1",
    fasta_headers = "p1",
    gene_symbols_or_id = "p1"
  )

  dataset = list(peptides = peptides, proteins = proteins, samples = samples, acquisition_mode = "dda")
  msdap::setup_contrasts(dataset, contrast_list = list(c("A", "B")))
}

get_contrast_col = function(dataset) {
  paste0("intensity_", dataset$contrasts[[1]]$label)
}

get_retained_peptides = function(peptides, col_intensity) {
  sort(unique(peptides$peptide_id[is.finite(peptides[[col_intensity]])]))
}

testthat::test_that("explicit NA preserves filter_dataset defaults", {
  dataset_default = make_groups_min_pass_dataset()
  dataset_na = make_groups_min_pass_dataset()

  result_default = suppressWarnings(msdap::filter_dataset(
    dataset_default,
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = TRUE,
    by_contrast = TRUE
  ))
  result_na = suppressWarnings(msdap::filter_dataset(
    dataset_na,
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = TRUE,
    by_contrast = TRUE,
    groups_min_pass = NA
  ))

  contrast_col = get_contrast_col(result_default)
  testthat::expect_equal(result_default$peptides$intensity_all_group, result_na$peptides$intensity_all_group)
  testthat::expect_equal(result_default$peptides[[contrast_col]], result_na$peptides[[contrast_col]])
})

testthat::test_that("groups_min_pass relaxes all_group filtering", {
  result_strict = suppressWarnings(msdap::filter_dataset(
    make_groups_min_pass_dataset(),
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = TRUE,
    by_contrast = FALSE
  ))
  result_relaxed = suppressWarnings(msdap::filter_dataset(
    make_groups_min_pass_dataset(),
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = TRUE,
    by_contrast = FALSE,
    groups_min_pass = 2
  ))

  testthat::expect_equal(get_retained_peptides(result_strict$peptides, "intensity_all_group"), "pep_all")
  testthat::expect_equal(get_retained_peptides(result_relaxed$peptides, "intensity_all_group"), c("pep_ab", "pep_all"))
})

testthat::test_that("groups_min_pass relaxes by_contrast filtering", {
  result_strict = suppressWarnings(msdap::filter_dataset(
    make_groups_min_pass_dataset(),
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = FALSE,
    by_contrast = TRUE
  ))
  result_relaxed = suppressWarnings(msdap::filter_dataset(
    make_groups_min_pass_dataset(),
    filter_min_detect = 2,
    norm_algorithm = "",
    by_group = FALSE,
    all_group = FALSE,
    by_contrast = TRUE,
    groups_min_pass = 1
  ))

  contrast_col = get_contrast_col(result_strict)
  testthat::expect_equal(get_retained_peptides(result_strict$peptides, contrast_col), c("pep_ab", "pep_all"))
  testthat::expect_equal(get_retained_peptides(result_relaxed$peptides, contrast_col), c("pep_a", "pep_ab", "pep_all"))
})

testthat::test_that("analysis_quickstart forwards groups_min_pass", {
  result = suppressWarnings(msdap::analysis_quickstart(
    make_groups_min_pass_dataset(),
    filter_min_detect = 2,
    filter_by_contrast = TRUE,
    norm_algorithm = "",
    dea_algorithm = character(0),
    diffdetect_min_peptides_observed = NA,
    diffdetect_min_samples_observed = NA,
    output_abundance_tables = FALSE,
    output_qc_report = FALSE,
    output_dir = NA,
    groups_min_pass = 1
  ))

  contrast_col = get_contrast_col(result)
  testthat::expect_true(contrast_col %in% colnames(result$peptides))
  testthat::expect_equal(get_retained_peptides(result$peptides, contrast_col), c("pep_a", "pep_ab", "pep_all"))
})
