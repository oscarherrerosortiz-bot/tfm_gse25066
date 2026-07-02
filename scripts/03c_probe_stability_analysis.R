#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(digest)
  library(yaml)
  library(glmnet)
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
annotation_path <- path("data", "processed", "GSE25066", "gpl96_probe_annotation.tsv.gz")
outer_path <- path("results", "predictive_protocol", "splits", "internal_outer_folds.tsv.gz")
inner_path <- path("results", "predictive_protocol", "splits", "internal_inner_folds.tsv.gz")
hyper_path <- path("results", "predictive_modeling", "tables", "phase3b_hyperparameter_selections.tsv")
phase3b_checks_path <- path("results", "predictive_modeling", "tables", "phase3b_checks.tsv")
config_path <- path("config", "predictive_protocol.yml")

out_dir <- path("results", "probe_stability")
tables_dir <- file.path(out_dir, "tables")
logs_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

console_log <- file.path(logs_dir, "phase3c_console.log")
if (file.exists(console_log)) invisible(file.remove(console_log))

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  cat(msg, "\n", file = console_log, append = TRUE)
}

expected_metadata_sha <- "0ac5539556d5856c9c43da9157a5c7af5fefa648e4603b35ba6c1bbe1f510bf3"
expected_expression_sha <- "fb2aafdc80bb8a136da4949a0ef9d78421628fb5701011e806dd7c95e1744a80"

status_priority <- c(
  "ANALYSIS_COMPLETE" = 0,
  "REVIEW_REQUIRED" = 1,
  "ANALYSIS_FAILED" = 2
)
final_status <- "ANALYSIS_COMPLETE"
warnings_seen <- character()

set_status <- function(status) {
  if (status_priority[[status]] > status_priority[[final_status]]) {
    final_status <<- status
  }
}

record_warning <- function(scope, warning_text) {
  warnings_seen <<- c(warnings_seen, sprintf("%s: %s", scope, warning_text))
}

sha256_file <- function(file) digest(file = file, algo = "sha256")

read_gz_tsv <- function(file) {
  as.data.table(read.delim(gzfile(file), sep = "\t", header = TRUE, check.names = FALSE))
}


available_memory_mb <- function() {
  meminfo <- "/proc/meminfo"
  if (file.exists(meminfo)) {
    lines <- readLines(meminfo, warn = FALSE)
    memtotal <- grep("^MemTotal:", lines, value = TRUE)
    if (length(memtotal) == 1L) {
      kb <- as.numeric(sub("^MemTotal:\\s+([0-9]+)\\s+kB.*$", "\\1", memtotal))
      if (is.finite(kb)) return(list(mb = kb / 1024, source = "proc_meminfo_memtotal"))
    }
  }
  free_out <- tryCatch(
    suppressWarnings(system2("free", "-m", stdout = TRUE, stderr = TRUE)),
    error = function(e) structure(character(), status = 127L)
  )
  free_status <- attr(free_out, "status")
  if ((is.null(free_status) || identical(free_status, 0L)) && length(free_out) >= 2L) {
    mem_line <- grep("^Mem:", free_out, value = TRUE)
    if (length(mem_line) == 1L) {
      parts <- strsplit(gsub("^Mem:\\s+", "", mem_line), "\\s+")[[1]]
      total_mb <- suppressWarnings(as.numeric(parts[[1]]))
      if (is.finite(total_mb)) return(list(mb = total_mb, source = "free_m_total"))
    }
  }
  out <- tryCatch(
    suppressWarnings(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = TRUE)),
    error = function(e) structure(character(), status = 127L)
  )
  sysctl_status <- attr(out, "status")
  if ((is.null(sysctl_status) || identical(sysctl_status, 0L)) && length(out) == 1L && grepl("^[0-9]+$", out)) {
    return(list(mb = as.numeric(out) / 1024^2, source = "sysctl_hw.memsize"))
  }
  list(mb = NA_real_, source = "unavailable")
}

write_session_info <- function() {
  writeLines(capture.output(sessionInfo()), file.path(logs_dir, "session_info.txt"))
}

stability_category <- function(freq) {
  fifelse(freq >= 0.70, "high",
    fifelse(freq >= 0.50, "moderate",
      fifelse(freq >= 0.20, "low", "occasional")
    )
  )
}

read_inputs <- function() {
  config <- yaml::read_yaml(config_path)
  hashes <- list(
    metadata = sha256_file(metadata_path),
    expression = sha256_file(expression_path),
    annotation = sha256_file(annotation_path),
    hyper = sha256_file(hyper_path)
  )
  metadata <- fread(metadata_path, sep = "\t", colClasses = "character")
  expression <- read_gz_tsv(expression_path)
  annotation <- read_gz_tsv(annotation_path)
  outer <- read_gz_tsv(outer_path)
  inner <- read_gz_tsv(inner_path)
  hyper <- fread(hyper_path)
  phase3b_checks <- fread(phase3b_checks_path)

  if (!identical(hashes$metadata, expected_metadata_sha) || !identical(hashes$expression, expected_expression_sha)) {
    set_status("ANALYSIS_FAILED")
    stop("metadata or expression hash mismatch")
  }
  if (!all(phase3b_checks$passed)) {
    set_status("REVIEW_REQUIRED")
    stop("las comprobaciones previas no son todas correctas")
  }
  if (nrow(metadata) != 488L || nrow(expression) != 488L || ncol(expression) != 22284L) {
    set_status("ANALYSIS_FAILED")
    stop("input dimensions do not match expected values")
  }
  if (!identical(metadata$sample_id, expression$sample_id)) {
    set_status("ANALYSIS_FAILED")
    stop("metadata and expression sample order mismatch")
  }
  primary <- metadata[primary_modeling_eligible == "True"]
  if (nrow(primary) != 482L) {
    set_status("ANALYSIS_FAILED")
    stop("primary population is not 482")
  }
  response_counts <- table(primary$response_pcr_vs_rd)
  if (response_counts[["pCR"]] != 98L || response_counts[["RD"]] != 384L) {
    set_status("ANALYSIS_FAILED")
    stop("primary response counts are not 98/384")
  }
  probe_ids <- colnames(expression)[-1L]
  if (length(probe_ids) != 22283L || anyDuplicated(probe_ids)) {
    set_status("ANALYSIS_FAILED")
    stop("probe IDs are invalid")
  }
  expr_primary <- expression[sample_id %in% primary$sample_id]
  setkey(expr_primary, sample_id)
  setkey(primary, sample_id)
  expr_primary <- expr_primary[primary$sample_id]
  x <- as.matrix(expr_primary[, ..probe_ids])
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) {
    set_status("ANALYSIS_FAILED")
    stop("expression contains missing or non-finite values")
  }
  primary[, response_binary := as.integer(response_pcr_vs_rd == "pCR")]
  primary[, stratum := paste(response_pcr_vs_rd, source_cohort, sep = "__")]
  if (!setequal(annotation$probe_id, probe_ids) || anyDuplicated(annotation$probe_id)) {
    set_status("REVIEW_REQUIRED")
    stop("annotation probe IDs do not match expression probes")
  }
  annotation <- annotation[match(probe_ids, probe_id)]

  input_integrity <- data.table(
    check = c(
      "metadata_sha256", "expression_sha256", "annotation_sha256",
      "hyperparameter_sha256",
      "metadata_rows", "expression_rows", "expression_probes",
      "primary_rows", "primary_pcr", "primary_rd", "annotation_rows",
      "phase3b_checks_all_true"
    ),
    observed = c(
      hashes$metadata, hashes$expression, hashes$annotation,
      hashes$hyper,
      nrow(metadata), nrow(expression), ncol(expression) - 1L,
      nrow(primary), response_counts[["pCR"]], response_counts[["RD"]],
      nrow(annotation), all(phase3b_checks$passed)
    ),
    expected = c(
      expected_metadata_sha, expected_expression_sha, "recorded",
      "recorded", 488L, 488L, 22283L, 482L, 98L, 384L, 22283L, TRUE
    ),
    passed = c(
      hashes$metadata == expected_metadata_sha,
      hashes$expression == expected_expression_sha,
      TRUE, TRUE,
      nrow(metadata) == 488L, nrow(expression) == 488L,
      ncol(expression) - 1L == 22283L, nrow(primary) == 482L,
      response_counts[["pCR"]] == 98L, response_counts[["RD"]] == 384L,
      nrow(annotation) == 22283L, all(phase3b_checks$passed)
    )
  )
  fwrite(input_integrity, file.path(tables_dir, "phase3c_input_integrity.tsv"), sep = "\t")

  list(
    config = config, metadata = metadata, primary = primary, x = x,
    probe_ids = probe_ids, annotation = annotation, outer = outer,
    inner = inner, hyper = hyper, hashes = hashes
  )
}

check_splits_and_hyper <- function(primary, outer, inner, hyper) {
  rows <- list()
  add <- function(scope, repeat_id, outer_fold, check, observed, expected, passed) {
    rows[[length(rows) + 1L]] <<- data.table(
      check_scope = scope,
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      check = check,
      observed = as.character(observed),
      expected = as.character(expected),
      passed = as.logical(passed)
    )
  }
  expected_folds <- CJ(repeat_id = 1:5, outer_fold = 1:5)
  trans_hyper <- hyper[validation_scheme == "internal" & model == "transcriptomic"]
  add("hyperparameter_file", NA_integer_, NA_integer_, "transcriptomic_internal_rows", nrow(trans_hyper), 25L, nrow(trans_hyper) == 25L)
  add("hyperparameter_file", NA_integer_, NA_integer_, "no_missing_alpha_lambda", all(is.finite(trans_hyper$selected_alpha) & is.finite(trans_hyper$selected_lambda)), TRUE, all(is.finite(trans_hyper$selected_alpha) & is.finite(trans_hyper$selected_lambda)))
  add("hyperparameter_file", NA_integer_, NA_integer_, "all_expected_repeat_outer_pairs", nrow(merge(expected_folds, trans_hyper[, .(repeat_id, outer_fold)], by = c("repeat_id", "outer_fold"))), 25L, nrow(merge(expected_folds, trans_hyper[, .(repeat_id, outer_fold)], by = c("repeat_id", "outer_fold"))) == 25L)

  no_outer_test_in_training <- TRUE
  inner_not_used_for_selection <- TRUE
  all_train_sizes_ok <- TRUE
  all_hyper_present <- TRUE
  for (i in seq_len(nrow(expected_folds))) {
    current_repeat <- expected_folds$repeat_id[[i]]
    current_outer <- expected_folds$outer_fold[[i]]
    test_ids <- outer[repeat_id == current_repeat & outer_fold == current_outer, sample_id]
    train_ids <- setdiff(primary$sample_id, test_ids)
    sub_inner <- inner[repeat_id == current_repeat & outer_fold == current_outer]
    hp <- trans_hyper[repeat_id == current_repeat & outer_fold == current_outer]
    no_outer_test_in_training <- no_outer_test_in_training && length(intersect(test_ids, train_ids)) == 0L
    inner_not_used_for_selection <- inner_not_used_for_selection && length(intersect(test_ids, sub_inner$sample_id)) == 0L
    all_train_sizes_ok <- all_train_sizes_ok && length(train_ids) + length(test_ids) == nrow(primary)
    all_hyper_present <- all_hyper_present && nrow(hp) == 1L
    add("outer_training", current_repeat, current_outer, "outer_test_not_used_for_probe_selection", length(intersect(test_ids, train_ids)), 0L, length(intersect(test_ids, train_ids)) == 0L)
    add("outer_training", current_repeat, current_outer, "outer_test_not_in_inner_file", length(intersect(test_ids, sub_inner$sample_id)), 0L, length(intersect(test_ids, sub_inner$sample_id)) == 0L)
    add("outer_training", current_repeat, current_outer, "hyperparameter_row_present", nrow(hp), 1L, nrow(hp) == 1L)
  }
  check_table <- rbindlist(rows)
  fwrite(check_table, file.path(tables_dir, "phase3c_split_usage_checks.tsv"), sep = "\t")
  list(
    check_table = check_table,
    no_outer_test_in_training = no_outer_test_in_training,
    inner_not_used_for_selection = inner_not_used_for_selection,
    all_train_sizes_ok = all_train_sizes_ok,
    all_hyper_present = all_hyper_present
  )
}

fit_fold_coefficients <- function(primary, x, probe_ids, outer, hyper, config) {
  trans_hyper <- hyper[validation_scheme == "internal" & model == "transcriptomic"]
  fold_rows <- list()
  hyper_checks <- list()
  runtime_rows <- list()
  warnings_local <- character()
  start_all <- proc.time()[["elapsed"]]

  for (current_repeat in sort(unique(trans_hyper$repeat_id))) {
    for (current_outer in sort(unique(trans_hyper$outer_fold))) {
      hp <- trans_hyper[repeat_id == current_repeat & outer_fold == current_outer]
      if (nrow(hp) != 1L) {
        set_status("REVIEW_REQUIRED")
        stop(sprintf("missing hyperparameter row for repeat %s outer %s", current_repeat, current_outer))
      }
      test_ids <- outer[repeat_id == current_repeat & outer_fold == current_outer, sample_id]
      train_idx <- which(!primary$sample_id %in% test_ids)
      y_train <- primary$response_binary[train_idx]
      x_train <- x[train_idx, , drop = FALSE]
      penalty <- rep(1, ncol(x_train))
      fold_start <- proc.time()[["elapsed"]]
      warn <- character()
      log_msg("FIT_START repeat=%s outer_fold=%s alpha=%s lambda=%s", current_repeat, current_outer, hp$selected_alpha, hp$selected_lambda)
      fit <- withCallingHandlers(
        glmnet(
          x = x_train,
          y = y_train,
          family = "binomial",
          alpha = hp$selected_alpha,
          lambda = hp$selected_lambda,
          penalty.factor = penalty,
          standardize = TRUE,
          intercept = TRUE,
          maxit = config$glmnet$maxit,
          thresh = config$glmnet$thresh
        ),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      elapsed <- proc.time()[["elapsed"]] - fold_start
      if (length(warn)) {
        warnings_local <- c(warnings_local, sprintf("repeat=%s outer=%s: %s", current_repeat, current_outer, paste(unique(warn), collapse = " | ")))
      }
      cf <- as.matrix(coef(fit, s = hp$selected_lambda))
      probe_coef <- cf[probe_ids, 1L]
      if (any(!is.finite(probe_coef))) {
        set_status("REVIEW_REQUIRED")
        stop(sprintf("non-finite coefficient repeat=%s outer=%s", current_repeat, current_outer))
      }
      nz <- which(abs(probe_coef) > 0)
      if (length(nz)) {
        fold_rows[[length(fold_rows) + 1L]] <- data.table(
          probe_id = probe_ids[nz],
          repeat_id = current_repeat,
          outer_fold = current_outer,
          coefficient = probe_coef[nz],
          coefficient_sign = fifelse(probe_coef[nz] > 0, "positive", "negative")
        )
      }
      hyper_checks[[length(hyper_checks) + 1L]] <- data.table(
        repeat_id = current_repeat,
        outer_fold = current_outer,
        alpha_reused = hp$selected_alpha,
        lambda_reused = hp$selected_lambda,
        n_penalized_nonzero_phase3b = hp$n_penalized_nonzero,
        n_nonzero_refit = length(nz),
        nonzero_count_difference = length(nz) - hp$n_penalized_nonzero,
        nonzero_count_matches_phase3b = length(nz) == hp$n_penalized_nonzero,
        outer_training_n = length(train_idx),
        outer_test_n = length(test_ids),
        all_sondes_penalized = all(penalty == 1),
        hyperparameters_changed = FALSE
      )
      runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
        model = "transcriptomic",
        repeat_id = current_repeat,
        outer_fold = current_outer,
        fit_time_sec = elapsed,
        memory_mb_approx = as.numeric(object.size(x_train)) / 1024^2,
        warnings = length(warn),
        nonfinite_coefficients = sum(!is.finite(probe_coef)),
        selected_probes = length(nz)
      )
      log_msg("FIT_END repeat=%s outer_fold=%s selected_probes=%s elapsed_sec=%.2f", current_repeat, current_outer, length(nz), elapsed)
    }
  }
  selected_long <- if (length(fold_rows)) rbindlist(fold_rows) else data.table(
    probe_id = character(), repeat_id = integer(), outer_fold = integer(),
    coefficient = numeric(), coefficient_sign = character()
  )
  hyper_checks_dt <- rbindlist(hyper_checks)
  runtime_dt <- rbindlist(runtime_rows)
  runtime_dt[, total_phase3c_fit_time_sec := proc.time()[["elapsed"]] - start_all]
  if (length(warnings_local)) record_warning("glmnet_refit", paste(unique(warnings_local), collapse = " | "))
  fwrite(hyper_checks_dt, file.path(tables_dir, "phase3c_hyperparameter_reuse_checks.tsv"), sep = "\t")
  list(selected_long = selected_long, hyper_checks = hyper_checks_dt, runtime = runtime_dt)
}

build_stability_table <- function(selected_long, annotation, probe_ids) {
  total_folds <- 25L
  if (nrow(selected_long)) {
    selected_long[, fold_record := sprintf("r%sf%s:%s:%.10g", repeat_id, outer_fold, coefficient_sign, coefficient)]
    agg <- selected_long[, .(
      selection_count = .N,
      positive_count = sum(coefficient_sign == "positive"),
      negative_count = sum(coefficient_sign == "negative"),
      mean_coefficient_selected_only = mean(coefficient),
      median_coefficient_selected_only = median(coefficient),
      mean_abs_coefficient_selected_only = mean(abs(coefficient)),
      max_abs_coefficient_selected_only = max(abs(coefficient)),
      selected_fold_records = paste(fold_record, collapse = ";")
    ), by = probe_id]
  } else {
    agg <- data.table(
      probe_id = character(), selection_count = integer(),
      positive_count = integer(), negative_count = integer(),
      mean_coefficient_selected_only = numeric(), median_coefficient_selected_only = numeric(),
      mean_abs_coefficient_selected_only = numeric(), max_abs_coefficient_selected_only = numeric(),
      selected_fold_records = character()
    )
  }
  stability <- data.table(probe_id = probe_ids)
  stability <- merge(stability, agg, by = "probe_id", all.x = TRUE, sort = FALSE)
  count_cols <- c("selection_count", "positive_count", "negative_count")
  for (cc in count_cols) stability[is.na(get(cc)), (cc) := 0L]
  stability[is.na(selected_fold_records), selected_fold_records := ""]
  stability[, selection_frequency := selection_count / total_folds]
  stability[, sign_consistency := fifelse(selection_count > 0, pmax(positive_count, negative_count) / selection_count, NA_real_)]
  stability[, stability_category := stability_category(selection_frequency)]
  stability[, interpretable_candidate := selection_frequency >= 0.50]
  stability <- merge(stability, annotation, by = "probe_id", all.x = TRUE, sort = FALSE)
  setcolorder(stability, c(
    "probe_id", "selection_count", "selection_frequency", "stability_category",
    "interpretable_candidate", "positive_count", "negative_count", "sign_consistency",
    "mean_coefficient_selected_only", "median_coefficient_selected_only",
    "mean_abs_coefficient_selected_only", "max_abs_coefficient_selected_only",
    "mapping_status", "gene_symbol_single", "gene_symbol_raw",
    "gene_symbols_normalized", "gene_title_raw", "gene_id_raw",
    setdiff(names(stability), c(
      "probe_id", "selection_count", "selection_frequency", "stability_category",
      "interpretable_candidate", "positive_count", "negative_count", "sign_consistency",
      "mean_coefficient_selected_only", "median_coefficient_selected_only",
      "mean_abs_coefficient_selected_only", "max_abs_coefficient_selected_only",
      "mapping_status", "gene_symbol_single", "gene_symbol_raw",
      "gene_symbols_normalized", "gene_title_raw", "gene_id_raw"
    ))
  ))
  stability[]
}

summaries_from_stability <- function(stability) {
  stability_summary <- stability[, .(
    n_probes = as.integer(.N),
    n_selected_at_least_once = as.integer(sum(selection_count > 0)),
    median_selection_count = as.numeric(median(selection_count)),
    max_selection_count = as.numeric(max(selection_count)),
    median_selection_frequency = median(selection_frequency),
    max_selection_frequency = max(selection_frequency),
    n_single_symbol = as.integer(sum(mapping_status == "single_symbol", na.rm = TRUE)),
    n_multiple_symbols = as.integer(sum(mapping_status == "multiple_symbols", na.rm = TRUE)),
    n_no_symbol = as.integer(sum(mapping_status == "no_symbol", na.rm = TRUE))
  ), by = stability_category]
  stability_summary[, category_order := match(stability_category, c("high", "moderate", "low", "occasional"))]
  setorder(stability_summary, category_order)
  stability_summary[, category_order := NULL]

  annotation_summary <- rbindlist(list(
    stability[, .(scope = "all_expression_probes", n_probes = .N), by = mapping_status],
    stability[selection_count > 0, .(scope = "selected_at_least_once", n_probes = .N), by = mapping_status],
    stability[selection_frequency >= 0.50, .(scope = "moderate_or_high_candidates", n_probes = .N), by = mapping_status]
  ), fill = TRUE)
  annotation_summary[is.na(mapping_status), mapping_status := "missing_annotation"]
  setcolorder(annotation_summary, c("scope", "mapping_status", "n_probes"))

  top_stable <- stability[selection_frequency >= 0.50]
  setorder(top_stable, -selection_frequency, -sign_consistency, -mean_abs_coefficient_selected_only, probe_id)
  top_stable <- top_stable[, .(
    probe_id, selection_count, selection_frequency, stability_category,
    positive_count, negative_count, sign_consistency,
    mean_coefficient_selected_only, median_coefficient_selected_only,
    mean_abs_coefficient_selected_only, max_abs_coefficient_selected_only,
    mapping_status, gene_symbol_single, gene_symbol_raw, gene_symbols_normalized,
    gene_title_raw, gene_id_raw
  )]
  list(stability_summary = stability_summary, annotation_summary = annotation_summary, top_stable = top_stable)
}

write_report <- function(final_status, inputs, split_checks, fit_out, stability, summaries, runtime) {
  selected_per_fold <- fit_out$hyper_checks[, .(
    min_selected = min(n_nonzero_refit),
    q1_selected = as.numeric(quantile(n_nonzero_refit, 0.25)),
    median_selected = median(n_nonzero_refit),
    q3_selected = as.numeric(quantile(n_nonzero_refit, 0.75)),
    max_selected = max(n_nonzero_refit),
    folds_without_selected_probes = sum(n_nonzero_refit == 0)
  )]
  high_n <- stability[stability_category == "high", .N]
  moderate_n <- stability[stability_category == "moderate", .N]
  low_n <- stability[stability_category == "low", .N]
  occasional_n <- stability[stability_category == "occasional", .N]
  candidate_text <- if ((high_n + moderate_n) > 0) {
    sprintf("Hay %s sondas con estabilidad moderada o alta candidatas a interpretacion prudente.", high_n + moderate_n)
  } else {
    "No hay sondas con estabilidad >= 0,50; no hay base suficiente para interpretacion robusta de sondas."
  }
  incident_lines <- character()
  if (length(warnings_seen)) incident_lines <- c(incident_lines, sprintf("- Warnings: %s", paste(unique(warnings_seen), collapse = " | ")))
  folds_zero <- fit_out$hyper_checks[n_nonzero_refit == 0]
  if (nrow(folds_zero)) {
    incident_lines <- c(incident_lines, sprintf("- Folds sin sondas seleccionadas: %s.", paste(sprintf("r%sf%s", folds_zero$repeat_id, folds_zero$outer_fold), collapse = ", ")))
  }
  count_mismatch <- fit_out$hyper_checks[nonzero_count_matches_phase3b == FALSE]
  if (nrow(count_mismatch)) {
    incident_lines <- c(incident_lines, sprintf(
      "- Recuentos no nulos distintos frente al diagnóstico previo en %s fold(s); alpha/lambda se reutilizaron sin cambios y la estabilidad se basa en el reajuste especificado en este análisis.",
      nrow(count_mismatch)
    ))
  }
  if (length(incident_lines) == 0L) incident_lines <- "- Sin incidencias relevantes registradas."

  report <- c(
    "# GSE25066: estabilidad de sondas",
    "",
    sprintf("Estado final: `%s`", final_status),
    "",
    "## Alcance",
    "",
    "Este análisis evalúa la estabilidad de seleccion de sondas del modelo transcriptomico en los 25 outer folds internos de la evaluación predictiva. Incluye anotacion GPL96 prudente. No incluye enriquecimiento funcional, GSEA/ORA, analisis diferencial ni interpretacion biologica fuerte.",
    "",
    "## Inputs",
    "",
    sprintf("- Metadatos: `%s`", rel(metadata_path)),
    sprintf("- Expresion: `%s`", rel(expression_path)),
    sprintf("- Anotacion GPL96: `%s`", rel(annotation_path)),
    sprintf("- Hiperparámetros reutilizados: `%s`", rel(hyper_path)),
    sprintf("- Hash metadatos: `%s`", inputs$hashes$metadata),
    sprintf("- Hash expresion: `%s`", inputs$hashes$expression),
    sprintf("- Hash anotacion: `%s`", inputs$hashes$annotation),
    sprintf("- Hash hiperparametros: `%s`", inputs$hashes$hyper),
    "- Dimensiones: 488 endpoint-known; 482 primarias; 98 pCR; 384 RD; 22.283 sondas.",
    "",
    "## Metodo",
    "",
    "- Se reutilizaron sin cambios los alpha/lambda seleccionados en la evaluación predictiva para el modelo transcriptomico.",
    "- Cada ajuste se rehizo solo sobre el outer training set correspondiente.",
    "- Se extrajeron coeficientes no nulos de sondas penalizadas.",
    "- La estabilidad se calculo como frecuencia de seleccion en 25 folds.",
    "- La anotacion GPL96 se unio por `probe_id`, preservando `single_symbol`, `multiple_symbols` y `no_symbol`.",
    "",
    "## Comprobaciones tecnicas",
    "",
    sprintf("- input integrity: %s", all(fread(file.path(tables_dir, "phase3c_input_integrity.tsv"))$passed)),
    sprintf("- split usage: %s", all(split_checks$check_table$passed)),
    sprintf("- hyperparameter reuse: %s", all(!fit_out$hyper_checks$hyperparameters_changed & is.finite(fit_out$hyper_checks$alpha_reused) & is.finite(fit_out$hyper_checks$lambda_reused))),
    "- no outer test used for probe selection: TRUE",
    "- population/probes unchanged: TRUE",
    "- annotation completed: TRUE",
    "- no DEG/GSEA/ORA: TRUE",
    "- no predictions by sample: TRUE",
    "",
    "## Resultados",
    "",
    "### Sondas seleccionadas por fold",
    "",
    paste(capture.output(print(selected_per_fold)), collapse = "\n"),
    "",
    "### Categorias de estabilidad",
    "",
    paste(capture.output(print(summaries$stability_summary)), collapse = "\n"),
    "",
    sprintf("- Alta estabilidad: %s sondas.", high_n),
    sprintf("- Estabilidad moderada: %s sondas.", moderate_n),
    sprintf("- Baja estabilidad: %s sondas.", low_n),
    sprintf("- Inestable/ocasional: %s sondas.", occasional_n),
    sprintf("- %s", candidate_text),
    "",
    "### Resumen de anotacion",
    "",
    paste(capture.output(print(summaries$annotation_summary)), collapse = "\n"),
    "",
    "## Limites de interpretacion",
    "",
    "Las sondas son la unidad principal. La anotacion a genes es auxiliar y conserva ambiguedades. Sondas de baja estabilidad u ocasionales no deben convertirse en historia biologica. No se realizo enriquecimiento funcional.",
    "",
    "## Incidencias",
    "",
    paste(incident_lines, collapse = "\n"),
    "",
    "## Archivos generados",
    "",
    "- `results/probe_stability/tables/phase3c_transcriptomic_probe_stability.tsv`",
    "- `results/probe_stability/tables/phase3c_transcriptomic_top_stable_probes.tsv`",
    "- `results/probe_stability/tables/phase3c_transcriptomic_annotation_summary.tsv`",
    "- `results/probe_stability/tables/phase3c_transcriptomic_stability_summary.tsv`",
    "- `results/probe_stability/tables/phase3c_hyperparameter_reuse_checks.tsv`",
    "- `results/probe_stability/logs/phase3c_console.log`",
    "",
    "## Recomendacion tecnica",
    "",
    if (final_status == "ANALYSIS_COMPLETE") "analysis_complete" else if (final_status == "REVIEW_REQUIRED") "review_required" else "analysis_failed"
  )
  writeLines(report, file.path(out_dir, "gse25066_phase3c_probe_stability_report.md"))
}


main <- function() {
  log_msg("STABILITY_START")
  inputs <- read_inputs()
  log_msg("input_hashes metadata=%s expression=%s annotation=%s hyper=%s", inputs$hashes$metadata, inputs$hashes$expression, inputs$hashes$annotation, inputs$hashes$hyper)
  log_msg("dimensions endpoint_known=488 primary=482 probes=22283")

  split_checks <- check_splits_and_hyper(inputs$primary, inputs$outer, inputs$inner, inputs$hyper)
  if (!all(split_checks$check_table$passed)) {
    set_status("ANALYSIS_FAILED")
    stop("split or hyperparameter check_table failed")
  }
  log_msg("split_usage_check=PASS")

  fit_out <- fit_fold_coefficients(inputs$primary, inputs$x, inputs$probe_ids, inputs$outer, inputs$hyper, inputs$config)
  if (any(fit_out$runtime$warnings > 0) || any(fit_out$runtime$nonfinite_coefficients > 0)) {
    set_status("REVIEW_REQUIRED")
  }

  stability <- build_stability_table(fit_out$selected_long, inputs$annotation, inputs$probe_ids)
  summaries <- summaries_from_stability(stability)
  fwrite(stability, file.path(tables_dir, "phase3c_transcriptomic_probe_stability.tsv"), sep = "\t")
  fwrite(summaries$top_stable, file.path(tables_dir, "phase3c_transcriptomic_top_stable_probes.tsv"), sep = "\t")
  fwrite(summaries$annotation_summary, file.path(tables_dir, "phase3c_transcriptomic_annotation_summary.tsv"), sep = "\t")
  fwrite(summaries$stability_summary, file.path(tables_dir, "phase3c_transcriptomic_stability_summary.tsv"), sep = "\t")

  memory_info <- available_memory_mb()
  memory_peak_mb <- max(fit_out$runtime$memory_mb_approx, na.rm = TRUE)
  memory_fraction <- if (is.finite(memory_info$mb)) memory_peak_mb / memory_info$mb else NA_real_
  memory_verified <- is.finite(memory_info$mb) && is.finite(memory_fraction) && memory_fraction < 0.70
  runtime <- copy(fit_out$runtime)
  runtime[, `:=`(
    available_memory_mb = memory_info$mb,
    memory_fraction_available = memory_fraction,
    memory_available_source = memory_info$source,
    memory_verification_passed = memory_verified
  )]
  runtime_summary <- runtime[, .(
    model = "summary",
    repeat_id = NA_integer_,
    outer_fold = NA_integer_,
    fit_time_sec = sum(fit_time_sec, na.rm = TRUE),
    memory_mb_approx = memory_peak_mb,
    warnings = sum(warnings, na.rm = TRUE),
    nonfinite_coefficients = sum(nonfinite_coefficients, na.rm = TRUE),
    selected_probes = sum(selected_probes, na.rm = TRUE),
    total_phase3c_fit_time_sec = max(total_phase3c_fit_time_sec, na.rm = TRUE),
    available_memory_mb = memory_info$mb,
    memory_fraction_available = memory_fraction,
    memory_available_source = memory_info$source,
    memory_verification_passed = memory_verified
  )]
  fwrite(rbindlist(list(runtime, runtime_summary), fill = TRUE), file.path(tables_dir, "phase3c_runtime_summary.tsv"), sep = "\t")
  if (!memory_verified) set_status("REVIEW_REQUIRED")

  high_n <- stability[stability_category == "high", .N]
  moderate_n <- stability[stability_category == "moderate", .N]
  low_n <- stability[stability_category == "low", .N]

  checks <- data.table(
    check_name = c(
      "input_integrity_passed",
      "split_usage_passed",
      "hyperparameter_reuse_passed",
      "no_outer_test_used_for_probe_selection",
      "population_unchanged",
      "probes_unchanged",
      "annotation_completed",
      "coefficients_finite",
      "memory_verification_passed",
      "no_deg_gsea_ora",
      "no_predictions_by_sample",
      "no_global_gene_collapse"
    ),
    passed = c(
      TRUE,
      all(split_checks$check_table$passed),
      all(!fit_out$hyper_checks$hyperparameters_changed) &&
        all(is.finite(fit_out$hyper_checks$alpha_reused)) &&
        all(is.finite(fit_out$hyper_checks$lambda_reused)) &&
        nrow(fit_out$hyper_checks) == 25L,
      TRUE,
      TRUE,
      nrow(stability) == 22283L,
      all(!is.na(stability$mapping_status)),
      all(fit_out$runtime$nonfinite_coefficients == 0),
      memory_verified,
      TRUE,
      TRUE,
      TRUE
    )
  )
  fwrite(checks, file.path(tables_dir, "phase3c_checks.tsv"), sep = "\t")
  if (!all(checks$passed)) set_status("REVIEW_REQUIRED")

  write_report(final_status, inputs, split_checks, fit_out, stability, summaries, runtime)
  write_session_info()

  log_msg("FINAL_STATUS=%s", final_status)
  log_msg("stability_counts high=%s moderate=%s low=%s occasional=%s", high_n, moderate_n, low_n, stability[stability_category == "occasional", .N])
  log_msg("files_written=%s", paste(c(
    rel(file.path(out_dir, "gse25066_phase3c_probe_stability_report.md")),
    rel(file.path(tables_dir, "phase3c_checks.tsv")),
    rel(file.path(tables_dir, "phase3c_transcriptomic_probe_stability.tsv")),
    rel(file.path(tables_dir, "phase3c_transcriptomic_top_stable_probes.tsv")),
    rel(console_log)
  ), collapse = ";"))
}

tryCatch(main(), error = function(e) {
  log_msg("ERROR=%s", conditionMessage(e))
  if (final_status == "ANALYSIS_COMPLETE") final_status <<- "REVIEW_REQUIRED"
  log_msg("FINAL_STATUS=%s", final_status)
  quit(status = 1L)
})
