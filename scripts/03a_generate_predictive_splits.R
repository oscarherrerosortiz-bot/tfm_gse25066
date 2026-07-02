#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(digest)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
ROOT <- if (!is.na(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(ROOT, "config", "predictive_protocol.yml"))) {
  ROOT <- normalizePath(getwd(), mustWork = TRUE)
}

path <- function(...) file.path(ROOT, ...)
rel <- function(p) sub(paste0("^", ROOT, "/?"), "", normalizePath(p, mustWork = FALSE))

metadata_path <- path("data", "processed", "GSE25066", "gse25066_analysis_metadata.tsv")
expression_path <- path("data", "processed", "GSE25066", "gse25066_expression_probe_endpoint_known.tsv.gz")
config_path <- path("config", "predictive_protocol.yml")

out_dir <- path("results", "predictive_protocol")
splits_dir <- path("results", "predictive_protocol", "splits")
tables_dir <- path("results", "predictive_protocol", "tables")
dir.create(splits_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

expected_metadata_sha <- "0ac5539556d5856c9c43da9157a5c7af5fefa648e4603b35ba6c1bbe1f510bf3"
expected_expression_sha <- "fb2aafdc80bb8a136da4949a0ef9d78421628fb5701011e806dd7c95e1744a80"

log_line <- function(...) cat(sprintf(...), "\n")

sha256_file <- function(file) digest(file = file, algo = "sha256")

stop_status <- function(status, message) {
  stop(sprintf("%s: %s", status, message), call. = FALSE)
}

read_inputs <- function() {
  if (!file.exists(metadata_path)) stop_status("INPUT_INTEGRITY_ERROR", "metadata file missing")
  if (!file.exists(expression_path)) stop_status("INPUT_INTEGRITY_ERROR", "expression file missing")
  if (!file.exists(config_path)) stop_status("REVIEW_REQUIRED_IMPLEMENTATION", "config file missing")

  config <- yaml::read_yaml(config_path)
  metadata_sha <- sha256_file(metadata_path)
  expression_sha <- sha256_file(expression_path)
  if (!identical(metadata_sha, expected_metadata_sha)) {
    stop_status("INPUT_INTEGRITY_ERROR", sprintf("metadata SHA mismatch: %s", metadata_sha))
  }
  if (!identical(expression_sha, expected_expression_sha)) {
    stop_status("INPUT_INTEGRITY_ERROR", sprintf("expression SHA mismatch: %s", expression_sha))
  }

  metadata <- fread(metadata_path, sep = "\t", colClasses = "character")
  expression <- as.data.table(read.delim(gzfile(expression_path), sep = "\t", header = TRUE, check.names = FALSE))
  if (nrow(metadata) != 488L) stop_status("INPUT_INTEGRITY_ERROR", "metadata row count is not 488")
  if (nrow(expression) != 488L || ncol(expression) != 22284L) {
    stop_status("INPUT_INTEGRITY_ERROR", sprintf("expression dimensions are %s x %s", nrow(expression), ncol(expression)))
  }
  if (anyDuplicated(metadata$sample_id)) stop_status("INPUT_INTEGRITY_ERROR", "duplicated metadata sample_id")
  if (anyDuplicated(expression$sample_id)) stop_status("INPUT_INTEGRITY_ERROR", "duplicated expression sample_id")
  if (!identical(expression$sample_id, metadata$sample_id)) stop_status("INPUT_INTEGRITY_ERROR", "metadata/expression sample order mismatch")

  primary <- metadata[primary_modeling_eligible == "True"]
  if (nrow(primary) != 482L) stop_status("INPUT_INTEGRITY_ERROR", "primary population is not 482")
  counts <- table(primary$response_pcr_vs_rd)
  if (!identical(as.integer(counts[["pCR"]]), 98L) || !identical(as.integer(counts[["RD"]]), 384L)) {
    stop_status("INPUT_INTEGRITY_ERROR", "primary response counts are not 98/384")
  }
  if ((ncol(expression) - 1L) != 22283L) stop_status("INPUT_INTEGRITY_ERROR", "probe count is not 22283")

  input_integrity <- data.table(
    check = c(
      "metadata_sha256", "expression_sha256", "metadata_rows", "expression_rows",
      "expression_probes", "primary_rows", "primary_pcr", "primary_rd"
    ),
    observed = c(
      metadata_sha, expression_sha, nrow(metadata), nrow(expression), ncol(expression) - 1L,
      nrow(primary), counts[["pCR"]], counts[["RD"]]
    ),
    expected = c(
      expected_metadata_sha, expected_expression_sha, 488L, 488L, 22283L, 482L, 98L, 384L
    ),
    passed = c(
      metadata_sha == expected_metadata_sha, expression_sha == expected_expression_sha,
      nrow(metadata) == 488L, nrow(expression) == 488L, ncol(expression) - 1L == 22283L,
      nrow(primary) == 482L, counts[["pCR"]] == 98L, counts[["RD"]] == 384L
    )
  )
  fwrite(input_integrity, path("results", "predictive_protocol", "tables", "input_integrity.tsv"), sep = "\t")
  list(config = config, metadata = metadata, expression = expression, primary = copy(primary))
}

make_stratified_folds <- function(dt, k, seed, stratum_col = "stratum") {
  set.seed(seed)
  folds <- integer(nrow(dt))
  for (stratum in sort(unique(dt[[stratum_col]]))) {
    idx <- which(dt[[stratum_col]] == stratum)
    if (length(idx) < k) {
      stop_status("REVIEW_REQUIRED_IMPLEMENTATION", sprintf("stratum %s has fewer than %s samples", stratum, k))
    }
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

validate_outer <- function(splits, primary) {
  balance_rows <- list()
  for (current_repeat_id in sort(unique(splits$repeat_id))) {
    sub <- splits[repeat_id == current_repeat_id]
    test_count <- sub[, .N, by = sample_id]
    if (!all(test_count$N == 1L)) stop_status("SPLIT_INTEGRITY_ERROR", "sample not exactly once as test per repeat")
    for (fold in sort(unique(sub$outer_fold))) {
      fold_dt <- sub[outer_fold == fold]
      if (anyDuplicated(fold_dt$sample_id)) stop_status("SPLIT_INTEGRITY_ERROR", "duplicated sample in outer fold")
      if (!all(c("pCR", "RD") %in% fold_dt$response)) {
        stop_status("REVIEW_REQUIRED_IMPLEMENTATION", "outer fold lacks pCR or RD")
      }
      if (length(unique(fold_dt$stratum)) != 8L) {
        stop_status("REVIEW_REQUIRED_IMPLEMENTATION", "outer fold lacks one of eight strata")
      }
      balance_row <- fold_dt[, .(n = .N), by = .(repeat_id, outer_fold, response, source_cohort, stratum)]
      balance_row[, split_type := "internal_outer"]
      balance_rows[[length(balance_rows) + 1L]] <- balance_row
    }
  }
  total_test <- splits[, .N, by = sample_id]
  if (!all(total_test$N == 5L)) stop_status("SPLIT_INTEGRITY_ERROR", "sample not exactly five times as test across repeats")
  rbindlist(balance_rows)
}

validate_inner <- function(inner_splits, outer_splits, primary) {
  balance_rows <- list()
  for (current_repeat_id in sort(unique(inner_splits$repeat_id))) {
    current_outer_folds <- sort(unique(inner_splits[repeat_id == current_repeat_id]$outer_fold))
    for (current_outer_fold in current_outer_folds) {
      outer_test <- outer_splits[repeat_id == current_repeat_id & outer_fold == current_outer_fold, sample_id]
      sub <- inner_splits[repeat_id == current_repeat_id & outer_fold == current_outer_fold]
      if (length(intersect(outer_test, sub$sample_id)) > 0L) {
        stop_status("SPLIT_INTEGRITY_ERROR", "outer test sample appears in inner split")
      }
      count <- sub[, .N, by = sample_id]
      if (!all(count$N == 1L)) stop_status("SPLIT_INTEGRITY_ERROR", "inner training sample not assigned exactly once")
      for (current_inner_fold in sort(unique(sub$inner_fold))) {
        val <- sub[inner_fold == current_inner_fold]
        if (!all(c("pCR", "RD") %in% val$response)) {
          stop_status("REVIEW_REQUIRED_IMPLEMENTATION", "inner fold lacks pCR or RD")
        }
        balance_row <- val[, .(n = .N), by = .(repeat_id, outer_fold, inner_fold, response, source_cohort, stratum)]
        balance_row[, split_type := "internal_inner"]
        balance_rows[[length(balance_rows) + 1L]] <- balance_row
      }
    }
  }
  rbindlist(balance_rows, fill = TRUE)
}

generate_internal_splits <- function(primary, config) {
  primary <- copy(primary)
  primary[, stratum := paste(response_pcr_vs_rd, source_cohort, sep = "__")]
  outer_rows <- list()
  seeds <- unlist(config$outer_seeds)
  for (i in seq_along(seeds)) {
    folds <- make_stratified_folds(primary, config$outer_folds, seeds[[i]])
    outer_rows[[i]] <- primary[
      ,
      .(
        sample_id,
        repeat_id = i,
        outer_seed = seeds[[i]],
        outer_fold = folds,
        response = response_pcr_vs_rd,
        source_cohort,
        stratum
      )
    ]
  }
  outer <- rbindlist(outer_rows)
  outer_balance <- validate_outer(outer, primary)

  inner_rows <- list()
  counter <- 0L
  for (current_repeat_id in sort(unique(outer$repeat_id))) {
    seed <- unique(outer[repeat_id == current_repeat_id]$outer_seed)
    for (current_outer_fold in sort(unique(outer$outer_fold))) {
      test_ids <- outer[repeat_id == current_repeat_id & outer_fold == current_outer_fold, sample_id]
      train <- primary[!sample_id %in% test_ids]
      folds <- make_stratified_folds(train, config$inner_folds, seed + current_outer_fold * 1000L)
      counter <- counter + 1L
      inner_rows[[counter]] <- train[
        ,
        .(
          repeat_id = current_repeat_id,
          outer_fold = current_outer_fold,
          sample_id,
          inner_fold = folds,
          response = response_pcr_vs_rd,
          source_cohort,
          stratum
        )
      ]
    }
  }
  inner <- rbindlist(inner_rows)
  inner_balance <- validate_inner(inner, outer, primary)
  list(outer = outer, inner = inner, balance = rbindlist(list(outer_balance, inner_balance), fill = TRUE))
}

generate_loco_splits <- function(primary, config) {
  primary <- copy(primary)
  primary[, stratum := paste(response_pcr_vs_rd, source_cohort, sep = "__")]
  domains <- primary[, .(sample_id, source_cohort, loco_test_domain = source_cohort, response = response_pcr_vs_rd, stratum)]
  inner_rows <- list()
  balance_rows <- list()
  counter <- 0L
  for (domain in sort(unique(primary$source_cohort))) {
    train <- primary[source_cohort != domain]
    folds <- make_stratified_folds(train, config$inner_folds, config$smoke_seed + match(domain, sort(unique(primary$source_cohort))) * 100L)
    counter <- counter + 1L
    inner_rows[[counter]] <- train[
      ,
      .(
        loco_test_domain = domain,
        sample_id,
        inner_fold = folds,
        response = response_pcr_vs_rd,
        source_cohort,
        stratum
      )
    ]
    balance_rows[[counter]] <- inner_rows[[counter]][
      ,
      .(n = .N),
      by = .(loco_test_domain, inner_fold, response, source_cohort, stratum)
    ]
    balance_rows[[counter]][, split_type := "loco_inner"]
  }
  list(domains = domains, inner = rbindlist(inner_rows), balance = rbindlist(balance_rows, fill = TRUE))
}

main <- function() {
  log_line("SPLIT_GENERATION_START")
  inputs <- read_inputs()
  primary <- inputs$primary
  config <- inputs$config
  internal <- generate_internal_splits(primary, config)
  loco <- generate_loco_splits(primary, config)

  fwrite(internal$outer, file.path(splits_dir, "internal_outer_folds.tsv.gz"), sep = "\t")
  fwrite(internal$inner, file.path(splits_dir, "internal_inner_folds.tsv.gz"), sep = "\t")
  fwrite(loco$domains, file.path(splits_dir, "loco_outer_domains.tsv"), sep = "\t")
  fwrite(loco$inner, file.path(splits_dir, "loco_inner_folds.tsv.gz"), sep = "\t")
  fwrite(rbindlist(list(internal$balance, loco$balance), fill = TRUE), file.path(tables_dir, "fold_balance.tsv"), sep = "\t")

  log_line("outer_folds_validated=TRUE")
  log_line("inner_folds_validated=TRUE")
  log_line("loco_structure_validated=TRUE")
  log_line("SPLIT_GENERATION_COMPLETE")
}

tryCatch(main(), error = function(e) {
  message(conditionMessage(e))
  quit(status = 1L)
})
