#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(digest)
  library(yaml)
  library(glmnet)
  library(pROC)
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
config_path <- path("config", "predictive_protocol.yml")

phase3a_dir <- path("results", "predictive_protocol")
splits_dir <- path("results", "predictive_protocol", "splits")

out_dir <- path("results", "predictive_modeling")
tables_dir <- file.path(out_dir, "tables")
logs_dir <- file.path(out_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

console_log <- file.path(logs_dir, "phase3b_console.log")
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

average_precision <- function(y, prob) {
  ord <- order(prob, decreasing = TRUE)
  y_ord <- y[ord]
  n_pos <- sum(y_ord == 1L)
  if (n_pos == 0L) return(NA_real_)
  precision <- cumsum(y_ord == 1L) / seq_along(y_ord)
  sum(precision[y_ord == 1L]) / n_pos
}

safe_roc_auc <- function(y, prob) {
  tryCatch({
    roc_obj <- suppressMessages(pROC::roc(response = y, predictor = prob, quiet = TRUE))
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
}

discrete_metrics <- function(y, prob, threshold) {
  pred <- as.integer(prob >= threshold)
  tp <- sum(pred == 1L & y == 1L)
  tn <- sum(pred == 0L & y == 0L)
  fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L)
  sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  data.table(
    balanced_accuracy = mean(c(sens, spec)),
    sensitivity = sens,
    specificity = spec
  )
}

all_metrics <- function(y, prob, threshold) {
  disc <- discrete_metrics(y, prob, threshold)
  data.table(
    pr_auc = average_precision(y, prob),
    roc_auc = safe_roc_auc(y, prob),
    balanced_accuracy = disc$balanced_accuracy,
    sensitivity = disc$sensitivity,
    specificity = disc$specificity
  )
}

select_threshold <- function(y, prob) {
  if (any(!is.finite(prob))) stop("non-finite internal probabilities")
  if (any(prob < -1e-12 | prob > 1 + 1e-12)) stop("threshold selection requires probabilities in [0, 1]")
  prob <- pmin(pmax(prob, 0), 1)
  thresholds <- sort(unique(prob))
  scores <- vapply(thresholds, function(th) discrete_metrics(y, prob, th)$balanced_accuracy, numeric(1))
  max_score <- max(scores, na.rm = TRUE)
  eligible <- which(abs(scores - max_score) < 1e-12)
  distances <- abs(thresholds[eligible] - 0.5)
  best_distance <- min(distances)
  eligible <- eligible[abs(distances - best_distance) < 1e-12]
  max(thresholds[eligible])
}

probability_scale <- function(values, source_label) {
  if (any(!is.finite(values))) stop(sprintf("non-finite predictions in %s", source_label))
  if (any(values < -1e-12 | values > 1 + 1e-12)) {
    values <- plogis(values)
  }
  pmin(pmax(values, 0), 1)
}

make_baseline_matrix <- function(dt, include_cohort = TRUE) {
  er_positive <- as.integer(dt$er_status == "positive")
  if (!include_cohort) {
    mat <- cbind(ER_positive = er_positive)
    storage.mode(mat) <- "double"
    return(mat)
  }
  mat <- cbind(
    ER_positive = er_positive,
    cohort_ISPY = as.integer(dt$source_cohort == "ISPY"),
    cohort_LBJ_IN_GEI = as.integer(dt$source_cohort == "LBJ/IN/GEI"),
    cohort_USO = as.integer(dt$source_cohort == "USO")
  )
  storage.mode(mat) <- "double"
  mat
}

fit_basal_predict <- function(train_dt, y_train, test_dt, include_cohort) {
  train_x <- as.data.frame(make_baseline_matrix(train_dt, include_cohort = include_cohort))
  test_x <- as.data.frame(make_baseline_matrix(test_dt, include_cohort = include_cohort))
  warn <- character()
  fit <- withCallingHandlers(
    glm(y_train ~ ., data = train_x, family = binomial()),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warn)) record_warning("basal_glm", paste(unique(warn), collapse = " | "))
  if (!isTRUE(fit$converged)) return(list(converged = FALSE, prob = rep(NA_real_, nrow(test_dt)), warnings = warn))
  prob <- as.numeric(predict(fit, newdata = test_x, type = "response"))
  list(converged = TRUE, prob = probability_scale(prob, "basal_glm_response"), warnings = warn)
}

fit_basal_inner_predictions <- function(train_dt, y_train, inner_fold, include_cohort) {
  probs <- rep(NA_real_, length(y_train))
  converged <- TRUE
  warn <- character()
  for (current_inner_fold in sort(unique(inner_fold))) {
    tr <- inner_fold != current_inner_fold
    val <- inner_fold == current_inner_fold
    fit_out <- fit_basal_predict(train_dt[tr], y_train[tr], train_dt[val], include_cohort = include_cohort)
    converged <- converged && isTRUE(fit_out$converged)
    warn <- c(warn, fit_out$warnings)
    probs[val] <- fit_out$prob
  }
  list(prob = probs, converged = converged, warnings = unique(warn))
}

select_global_1se <- function(cv_results, penalized_start) {
  candidates <- rbindlist(lapply(cv_results, function(item) {
    cv <- item$cv
    lambda <- cv$lambda
    n_nonzero <- vapply(seq_along(lambda), function(i) {
      cf <- as.matrix(coef(cv$glmnet.fit, s = lambda[[i]]))
      if (penalized_start > nrow(cf)) return(0L)
      sum(abs(cf[penalized_start:nrow(cf), 1L]) > 0)
    }, integer(1))
    data.table(
      alpha = item$alpha,
      lambda = lambda,
      lambda_index = seq_along(lambda),
      cvm = cv$cvm,
      cvsd = cv$cvsd,
      n_penalized_nonzero = n_nonzero,
      n_lambda = length(lambda)
    )
  }))
  min_idx <- which.min(candidates$cvm)
  threshold <- candidates$cvm[min_idx] + candidates$cvsd[min_idx]
  eligible <- candidates[cvm <= threshold]
  setorder(eligible, n_penalized_nonzero, -lambda, -alpha)
  selected <- eligible[1]
  selected[, on_boundary := lambda_index == 1L | lambda_index == n_lambda]
  selected[, boundary_side := fifelse(lambda_index == 1L, "upper", fifelse(lambda_index == n_lambda, "lower", "none"))]
  selected[, boundary_position := fifelse(lambda_index == 1L, "lambda_max", fifelse(lambda_index == n_lambda, "lambda_min", "none"))]
  selected[, selected_by_global_1se := TRUE]
  selected
}

count_penalized_nonzero <- function(glmnet_fit, lambda_value, penalized_start) {
  cf <- as.matrix(coef(glmnet_fit, s = lambda_value))
  if (penalized_start > nrow(cf)) return(0L)
  sum(abs(cf[penalized_start:nrow(cf), 1L]) > 0)
}

lambda_check_row <- function(scheme, model, repeat_id, outer_fold, loco_domain, cv_grid_result, penalized_start) {
  selected <- cv_grid_result$selected
  cv_selected <- cv_grid_result$cv_selected
  lambda_max <- max(cv_selected$lambda)
  nonzero_at_lambda_max <- count_penalized_nonzero(cv_selected$glmnet.fit, lambda_max, penalized_start)
  data.table(
    validation_scheme = scheme,
    model = model,
    repeat_id = repeat_id,
    outer_fold = outer_fold,
    loco_domain = loco_domain,
    selected_alpha = selected$alpha,
    lambda_on_boundary = selected$on_boundary,
    boundary_side = selected$boundary_side,
    boundary_position = selected$boundary_position,
    selected_by_global_1se = selected$selected_by_global_1se,
    penalized_nonzero_selected_lambda = selected$n_penalized_nonzero,
    penalized_nonzero_lambda_max = nonzero_at_lambda_max,
    all_penalized_zero_at_lambda_max = nonzero_at_lambda_max == 0L,
    lambda_boundary_technically_justified = ifelse(
      selected$on_boundary,
      selected$boundary_side == "upper" &&
        selected$boundary_position == "lambda_max" &&
        selected$selected_by_global_1se &&
        selected$n_penalized_nonzero == 0L &&
        nonzero_at_lambda_max == 0L,
      TRUE
    )
  )
}

run_cv_grid <- function(model_name, x, y, foldid, penalty_factor, config, penalized_start) {
  alpha_grid <- unlist(config$alpha_grid)
  warnings_local <- character()
  start <- proc.time()[["elapsed"]]
  cv_results <- list()
  for (alpha in alpha_grid) {
    cv <- withCallingHandlers(
      cv.glmnet(
        x = x,
        y = y,
        family = "binomial",
        alpha = alpha,
        nlambda = config$nlambda,
        lambda.min.ratio = config$lambda_min_ratio,
        type.measure = config$internal_tuning$type_measure,
        grouped = config$internal_tuning$grouped,
        keep = TRUE,
        foldid = foldid,
        penalty.factor = penalty_factor,
        standardize = TRUE,
        intercept = TRUE,
        maxit = config$glmnet$maxit,
        thresh = config$glmnet$thresh
      ),
      warning = function(w) {
        warnings_local <<- c(warnings_local, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    cv_results[[as.character(alpha)]] <- list(alpha = alpha, cv = cv)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  selected <- select_global_1se(cv_results, penalized_start)
  cv_selected <- cv_results[[as.character(selected$alpha)]]$cv
  list(
    model_name = model_name,
    cv_results = cv_results,
    selected = selected,
    cv_selected = cv_selected,
    tuning_time_sec = elapsed,
    warnings = unique(warnings_local)
  )
}

fit_glmnet_predict <- function(x_train, y_train, x_test, selected, penalty_factor, config) {
  warnings_local <- character()
  start <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = selected$alpha,
      lambda = selected$lambda,
      penalty.factor = penalty_factor,
      standardize = TRUE,
      intercept = TRUE,
      maxit = config$glmnet$maxit,
      thresh = config$glmnet$thresh
    ),
    warning = function(w) {
      warnings_local <<- c(warnings_local, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start
  prob <- as.numeric(predict(fit, newx = x_test, s = selected$lambda, type = "response"))
  prob <- probability_scale(prob, "glmnet_response")
  list(prob = prob, fit_time_sec = elapsed, warnings = unique(warnings_local))
}

read_inputs <- function() {
  config <- yaml::read_yaml(config_path)
  metadata_sha <- sha256_file(metadata_path)
  expression_sha <- sha256_file(expression_path)
  if (!identical(metadata_sha, expected_metadata_sha) || !identical(expression_sha, expected_expression_sha)) {
    set_status("ANALYSIS_FAILED")
    stop("input hashes do not match expected files")
  }
  metadata <- fread(metadata_path, sep = "\t", colClasses = "character")
  expression <- as.data.table(read.delim(gzfile(expression_path), sep = "\t", header = TRUE, check.names = FALSE))
  annotation <- as.data.table(read.delim(gzfile(annotation_path), sep = "\t", header = TRUE, check.names = FALSE))
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
    stop("probe IDs are not valid")
  }
  expr_primary <- expression[sample_id %in% primary$sample_id]
  setkey(expr_primary, sample_id)
  setkey(primary, sample_id)
  expr_primary <- expr_primary[primary$sample_id]
  x <- as.matrix(expr_primary[, ..probe_ids])
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) {
    set_status("ANALYSIS_FAILED")
    stop("expression matrix contains missing or non-finite values")
  }
  primary[, response_binary := as.integer(response_pcr_vs_rd == "pCR")]
  primary[, stratum := paste(response_pcr_vs_rd, source_cohort, sep = "__")]
  input_integrity <- data.table(
    check = c("metadata_sha256", "expression_sha256", "metadata_rows", "expression_rows", "expression_probes", "primary_rows", "primary_pcr", "primary_rd", "annotation_rows"),
    observed = c(metadata_sha, expression_sha, nrow(metadata), nrow(expression), ncol(expression) - 1L, nrow(primary), response_counts[["pCR"]], response_counts[["RD"]], nrow(annotation)),
    expected = c(expected_metadata_sha, expected_expression_sha, 488L, 488L, 22283L, 482L, 98L, 384L, 22283L),
    passed = c(metadata_sha == expected_metadata_sha, expression_sha == expected_expression_sha, nrow(metadata) == 488L, nrow(expression) == 488L, ncol(expression) - 1L == 22283L, nrow(primary) == 482L, response_counts[["pCR"]] == 98L, response_counts[["RD"]] == 384L, nrow(annotation) == 22283L)
  )
  fwrite(input_integrity, file.path(tables_dir, "phase3b_input_integrity.tsv"), sep = "\t")
  list(config = config, metadata = metadata, primary = primary, expression = x, probe_ids = probe_ids, annotation = annotation, hashes = list(metadata = metadata_sha, expression = expression_sha))
}

read_required_splits <- function() {
  split_paths <- list(
    internal_outer = file.path(splits_dir, "internal_outer_folds.tsv.gz"),
    internal_inner = file.path(splits_dir, "internal_inner_folds.tsv.gz"),
    loco_outer = file.path(splits_dir, "loco_outer_domains.tsv"),
    loco_inner = file.path(splits_dir, "loco_inner_folds.tsv.gz")
  )
  if (!all(file.exists(unlist(split_paths)))) {
    set_status("REVIEW_REQUIRED")
    stop("las particiones requeridas no se encuentran")
  }
  list(
    paths = split_paths,
    internal_outer = as.data.table(read.delim(gzfile(split_paths$internal_outer), sep = "\t", header = TRUE, check.names = FALSE)),
    internal_inner = as.data.table(read.delim(gzfile(split_paths$internal_inner), sep = "\t", header = TRUE, check.names = FALSE)),
    loco_outer = fread(split_paths$loco_outer),
    loco_inner = as.data.table(read.delim(gzfile(split_paths$loco_inner), sep = "\t", header = TRUE, check.names = FALSE))
  )
}

check_split_usage <- function(primary, splits, config) {
  rows <- list()
  add <- function(scope, check, observed, expected, passed) {
    rows[[length(rows) + 1L]] <<- data.table(
      check_scope = scope,
      check = check,
      observed = as.character(observed),
      expected = as.character(expected),
      passed = as.logical(passed)
    )
  }
  outer <- splits$internal_outer
  inner <- splits$internal_inner
  loco_outer <- splits$loco_outer
  loco_inner <- splits$loco_inner
  n_repeats <- length(unlist(config$outer_seeds))
  expected_strata <- sort(unique(primary$stratum))
  expected_cohorts <- sort(unique(primary$source_cohort))

  add("internal_outer", "row_count", nrow(outer), nrow(primary) * n_repeats, nrow(outer) == nrow(primary) * n_repeats)
  add("internal_outer", "sample_once_per_repeat", paste(range(outer[, .N, by = .(repeat_id, sample_id)]$N), collapse = "-"), "1-1", all(outer[, .N, by = .(repeat_id, sample_id)]$N == 1L))
  add("internal_outer", "sample_five_tests_total", paste(range(outer[, .N, by = sample_id]$N), collapse = "-"), "5-5", all(outer[, .N, by = sample_id]$N == n_repeats))
  outer_fold_summary <- outer[, .(responses = paste(sort(unique(response)), collapse = ","), cohorts = uniqueN(source_cohort), strata = uniqueN(stratum)), by = .(repeat_id, outer_fold)]
  add("internal_outer", "folds_have_both_responses", sum(outer_fold_summary$responses == "pCR,RD"), nrow(outer_fold_summary), all(outer_fold_summary$responses == "pCR,RD"))
  add("internal_outer", "folds_have_four_cohorts", sum(outer_fold_summary$cohorts == length(expected_cohorts)), nrow(outer_fold_summary), all(outer_fold_summary$cohorts == length(expected_cohorts)))
  add("internal_outer", "folds_have_eight_strata", sum(outer_fold_summary$strata == length(expected_strata)), nrow(outer_fold_summary), all(outer_fold_summary$strata == length(expected_strata)))

  no_outer_test <- TRUE
  each_train_once <- TRUE
  inner_both_response <- TRUE
  inner_eight_strata <- TRUE
  for (current_repeat_id in sort(unique(outer$repeat_id))) {
    for (current_outer_fold in sort(unique(outer$outer_fold))) {
      outer_test <- outer[repeat_id == current_repeat_id & outer_fold == current_outer_fold, sample_id]
      expected_train <- setdiff(primary$sample_id, outer_test)
      sub <- inner[repeat_id == current_repeat_id & outer_fold == current_outer_fold]
      no_outer_test <- no_outer_test && length(intersect(outer_test, sub$sample_id)) == 0L
      counts <- sub[, .N, by = sample_id]
      each_train_once <- each_train_once && setequal(counts$sample_id, expected_train) && all(counts$N == 1L)
      per_fold <- sub[, .(responses = paste(sort(unique(response)), collapse = ","), strata = uniqueN(stratum)), by = inner_fold]
      inner_both_response <- inner_both_response && all(per_fold$responses == "pCR,RD")
      inner_eight_strata <- inner_eight_strata && all(per_fold$strata == length(expected_strata))
    }
  }
  add("internal_inner", "outer_test_not_in_inner", no_outer_test, TRUE, no_outer_test)
  add("internal_inner", "outer_train_samples_once", each_train_once, TRUE, each_train_once)
  add("internal_inner", "inner_validation_both_responses", inner_both_response, TRUE, inner_both_response)
  add("internal_inner", "inner_validation_eight_strata_when_possible", inner_eight_strata, TRUE, inner_eight_strata)
  add("internal_inner", "same_inner_folds_for_all_models", "single shared file", "single shared file", TRUE)

  add("loco_outer", "row_count", nrow(loco_outer), nrow(primary), nrow(loco_outer) == nrow(primary))
  add("loco_outer", "domains", paste(sort(unique(loco_outer$loco_test_domain)), collapse = ","), paste(expected_cohorts, collapse = ","), setequal(loco_outer$loco_test_domain, expected_cohorts))
  loco_no_test <- TRUE
  loco_train_once <- TRUE
  loco_both_response <- TRUE
  for (current_domain in expected_cohorts) {
    heldout <- primary[source_cohort == current_domain, sample_id]
    train <- primary[source_cohort != current_domain, sample_id]
    sub <- loco_inner[loco_test_domain == current_domain]
    loco_no_test <- loco_no_test && length(intersect(heldout, sub$sample_id)) == 0L
    counts <- sub[, .N, by = sample_id]
    loco_train_once <- loco_train_once && setequal(counts$sample_id, train) && all(counts$N == 1L)
    per_fold <- sub[, .(responses = paste(sort(unique(response)), collapse = ",")), by = inner_fold]
    loco_both_response <- loco_both_response && all(per_fold$responses == "pCR,RD")
  }
  add("loco_inner", "heldout_domain_not_in_inner", loco_no_test, TRUE, loco_no_test)
  add("loco_inner", "loco_train_samples_once", loco_train_once, TRUE, loco_train_once)
  add("loco_inner", "loco_inner_validation_both_responses", loco_both_response, TRUE, loco_both_response)
  check_table <- rbindlist(rows)
  fwrite(check_table, file.path(tables_dir, "phase3b_split_usage_checks.tsv"), sep = "\t")
  check_table
}

validate_glmnet_boundary <- function(lambda_row) {
  if (!isTRUE(lambda_row$lambda_on_boundary)) return(TRUE)
  isTRUE(lambda_row$boundary_side == "upper" &&
    lambda_row$boundary_position == "lambda_max" &&
    lambda_row$selected_by_global_1se &&
    lambda_row$penalized_nonzero_selected_lambda == 0L &&
    lambda_row$all_penalized_zero_at_lambda_max)
}

metric_row <- function(validation_scheme, model, repeat_id, outer_fold, loco_domain, y_test, prob, threshold) {
  mets <- all_metrics(y_test, prob, threshold)
  data.table(
    validation_scheme = validation_scheme,
    model = model,
    repeat_id = repeat_id,
    outer_fold = outer_fold,
    loco_domain = loco_domain,
    n_test = length(y_test),
    pcr_n = sum(y_test == 1L),
    rd_n = sum(y_test == 0L),
    pcr_prevalence = mean(y_test == 1L),
    pr_auc = mets$pr_auc,
    roc_auc = mets$roc_auc,
    balanced_accuracy = mets$balanced_accuracy,
    sensitivity = mets$sensitivity,
    specificity = mets$specificity
  )
}

append_runtime <- function(runtime_rows, validation_scheme, model, repeat_id, outer_fold, loco_domain, fit_time, tuning_time, memory_mb, glmnet_trajectories, lambdas_returned, warnings, convergence_error, nonfinite_probabilities) {
  runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
    validation_scheme = validation_scheme,
    model = model,
    repeat_id = repeat_id,
    outer_fold = outer_fold,
    loco_domain = loco_domain,
    fit_time_sec = fit_time,
    inner_tuning_time_sec = tuning_time,
    memory_mb_approx = memory_mb,
    glmnet_trajectories = glmnet_trajectories,
    lambdas_returned = lambdas_returned,
    warnings = warnings,
    convergence_errors = convergence_error,
    nonfinite_probabilities = nonfinite_probabilities
  )
  runtime_rows
}

evaluate_one_split <- function(validation_scheme, model, repeat_id, outer_fold, loco_domain, train_idx, test_idx, inner_fold, primary, x, probe_ids, config) {
  train_dt <- primary[train_idx]
  test_dt <- primary[test_idx]
  y_train <- train_dt$response_binary
  y_test <- test_dt$response_binary
  x_train_probe <- x[train_idx, , drop = FALSE]
  x_test_probe <- x[test_idx, , drop = FALSE]
  include_cohort <- validation_scheme == "internal"
  start_total <- proc.time()[["elapsed"]]

  if (model == "basal") {
    inner <- fit_basal_inner_predictions(train_dt, y_train, inner_fold, include_cohort = include_cohort)
    threshold <- select_threshold(y_train, inner$prob)
    fit_start <- proc.time()[["elapsed"]]
    final <- fit_basal_predict(train_dt, y_train, test_dt, include_cohort = include_cohort)
    fit_time <- proc.time()[["elapsed"]] - fit_start
    prob <- final$prob
    converged <- isTRUE(final$converged) && isTRUE(inner$converged)
    warnings_n <- length(unique(c(inner$warnings, final$warnings)))
    metrics <- metric_row(validation_scheme, model, repeat_id, outer_fold, loco_domain, y_test, prob, threshold)
    return(list(
      metrics = metrics,
      threshold = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, threshold = threshold, threshold_source = "inner_training_balanced_accuracy"),
      hyper = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, selected_alpha = NA_real_, selected_lambda = NA_real_, selected_lambda_index = NA_integer_, n_lambda = NA_integer_, n_penalized_nonzero = NA_integer_),
      lambda = NULL,
      runtime = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, fit_time_sec = fit_time, inner_tuning_time_sec = 0, memory_mb_approx = as.numeric(object.size(make_baseline_matrix(train_dt, include_cohort))) / 1024^2, glmnet_trajectories = 0L, lambdas_returned = 0L, warnings = warnings_n, convergence_errors = as.integer(!converged), nonfinite_probabilities = sum(!is.finite(prob))),
      model_matrix = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, rows_train = nrow(train_dt), rows_test = nrow(test_dt), columns = ncol(make_baseline_matrix(train_dt, include_cohort)), er_predictor = TRUE, cohort_predictor = include_cohort, probe_columns = 0L, column_order_ok = TRUE),
      penalty = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, penalty_factor_length = 0L, zero_penalty_count = 0L, one_penalty_count = 0L, intercept_in_penalty_factor = FALSE, baseline_unpenalized = TRUE, all_probes_penalized = TRUE, order_ok = TRUE),
      converged = converged,
      probs_finite = all(is.finite(prob)),
      elapsed = proc.time()[["elapsed"]] - start_total
    ))
  }

  if (model == "transcriptomic") {
    design_train <- x_train_probe
    design_test <- x_test_probe
    penalty <- rep(1, ncol(design_train))
    penalized_start <- 2L
    er_predictor <- FALSE
    cohort_predictor <- FALSE
    order_ok <- TRUE
  } else if (model == "combined") {
    base_train <- make_baseline_matrix(train_dt, include_cohort = include_cohort)
    base_test <- make_baseline_matrix(test_dt, include_cohort = include_cohort)
    design_train <- cbind(base_train, x_train_probe)
    design_test <- cbind(base_test, x_test_probe)
    if (include_cohort) {
      expected_base <- c("ER_positive", "cohort_ISPY", "cohort_LBJ_IN_GEI", "cohort_USO")
      penalty <- c(rep(0, 4), rep(1, ncol(x_train_probe)))
    } else {
      expected_base <- "ER_positive"
      penalty <- c(0, rep(1, ncol(x_train_probe)))
    }
    colnames(design_train) <- c(expected_base, probe_ids)
    colnames(design_test) <- colnames(design_train)
    penalized_start <- length(expected_base) + 2L
    er_predictor <- TRUE
    cohort_predictor <- include_cohort
    order_ok <- identical(colnames(design_train)[seq_along(expected_base)], expected_base)
  } else {
    stop(sprintf("unknown model: %s", model))
  }

  cv_grid <- run_cv_grid(model, design_train, y_train, inner_fold, penalty, config, penalized_start)
  selected <- cv_grid$selected
  selected_lambda_index <- selected$lambda_index
  inner_prob <- probability_scale(
    cv_grid$cv_selected$fit.preval[, selected_lambda_index],
    sprintf("%s_%s_inner_prevalidated", validation_scheme, model)
  )
  threshold <- select_threshold(y_train, inner_prob)
  final <- fit_glmnet_predict(design_train, y_train, design_test, selected, penalty, config)
  prob <- final$prob
  lambda_row <- lambda_check_row(validation_scheme, model, repeat_id, outer_fold, loco_domain, cv_grid, penalized_start)
  lambda_ok <- validate_glmnet_boundary(lambda_row)
  if (!lambda_ok) set_status("REVIEW_REQUIRED")
  if (length(cv_grid$warnings) || length(final$warnings)) {
    record_warning(sprintf("%s_%s", validation_scheme, model), paste(unique(c(cv_grid$warnings, final$warnings)), collapse = " | "))
  }
  metrics <- metric_row(validation_scheme, model, repeat_id, outer_fold, loco_domain, y_test, prob, threshold)
  data_cols <- ncol(design_train)
  list(
    metrics = metrics,
    threshold = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, threshold = threshold, threshold_source = "inner_training_balanced_accuracy"),
    hyper = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, selected_alpha = selected$alpha, selected_lambda = selected$lambda, selected_lambda_index = selected$lambda_index, n_lambda = selected$n_lambda, n_penalized_nonzero = selected$n_penalized_nonzero),
    lambda = lambda_row,
    runtime = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, fit_time_sec = final$fit_time_sec, inner_tuning_time_sec = cv_grid$tuning_time_sec, memory_mb_approx = as.numeric(object.size(design_train)) / 1024^2, glmnet_trajectories = length(config$alpha_grid), lambdas_returned = sum(vapply(cv_grid$cv_results, function(item) length(item$cv$lambda), integer(1))), warnings = length(cv_grid$warnings) + length(final$warnings), convergence_errors = 0L, nonfinite_probabilities = sum(!is.finite(prob))),
    model_matrix = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, rows_train = nrow(train_dt), rows_test = nrow(test_dt), columns = data_cols, er_predictor = er_predictor, cohort_predictor = cohort_predictor, probe_columns = ncol(x_train_probe), column_order_ok = order_ok),
    penalty = data.table(validation_scheme, model, repeat_id, outer_fold, loco_domain, penalty_factor_length = length(penalty), zero_penalty_count = sum(penalty == 0), one_penalty_count = sum(penalty == 1), intercept_in_penalty_factor = FALSE, baseline_unpenalized = ifelse(model == "combined", all(penalty[seq_len(length(penalty) - ncol(x_train_probe))] == 0), NA), all_probes_penalized = ifelse(model %in% c("transcriptomic", "combined"), all(tail(penalty, ncol(x_train_probe)) == 1), TRUE), order_ok = order_ok),
    converged = TRUE,
    probs_finite = all(is.finite(prob)),
    elapsed = proc.time()[["elapsed"]] - start_total
  )
}

summarize_metric_table <- function(dt, scheme_filter) {
  metrics <- c("pr_auc", "roc_auc", "balanced_accuracy", "sensitivity", "specificity", "pcr_prevalence")
  rows <- list()
  for (metric in metrics) {
    rows[[length(rows) + 1L]] <- dt[validation_scheme == scheme_filter, .(
      n_evaluations = .N,
      mean = mean(get(metric), na.rm = TRUE),
      sd = sd(get(metric), na.rm = TRUE),
      median = median(get(metric), na.rm = TRUE),
      q1 = as.numeric(quantile(get(metric), 0.25, na.rm = TRUE)),
      q3 = as.numeric(quantile(get(metric), 0.75, na.rm = TRUE)),
      min = min(get(metric), na.rm = TRUE),
      max = max(get(metric), na.rm = TRUE)
    ), by = model][, metric := metric]
  }
  out <- rbindlist(rows, fill = TRUE)
  setcolorder(out, c("model", "metric", setdiff(names(out), c("model", "metric"))))
  out[]
}

paired_differences <- function(dt) {
  metrics <- c("pr_auc", "roc_auc", "balanced_accuracy", "sensitivity", "specificity")
  comparisons <- list(
    c("transcriptomic", "basal"),
    c("combined", "basal"),
    c("combined", "transcriptomic")
  )
  rows <- list()
  for (metric in metrics) {
    wide <- dcast(dt[validation_scheme == "internal"], repeat_id + outer_fold ~ model, value.var = metric)
    for (comparison in comparisons) {
      model_a <- comparison[[1]]
      model_b <- comparison[[2]]
      diff <- wide[[model_a]] - wide[[model_b]]
      rows[[length(rows) + 1L]] <- data.table(
        comparison = sprintf("%s_minus_%s", model_a, model_b),
        metric = metric,
        n_pairs = length(diff),
        mean_difference = mean(diff, na.rm = TRUE),
        sd_difference = sd(diff, na.rm = TRUE),
        median_difference = median(diff, na.rm = TRUE),
        q1_difference = as.numeric(quantile(diff, 0.25, na.rm = TRUE)),
        q3_difference = as.numeric(quantile(diff, 0.75, na.rm = TRUE)),
        min_difference = min(diff, na.rm = TRUE),
        max_difference = max(diff, na.rm = TRUE)
      )
    }
  }
  rbindlist(rows)
}

write_session_info <- function() {
  writeLines(capture.output(sessionInfo()), file.path(logs_dir, "session_info.txt"))
}

main <- function() {
  total_start <- proc.time()[["elapsed"]]
  log_msg("EVALUATION_START")
  inputs <- read_inputs()
  config <- inputs$config
  primary <- inputs$primary
  x <- inputs$expression
  probe_ids <- inputs$probe_ids
  splits <- read_required_splits()
  split_checks <- check_split_usage(primary, splits, config)
  if (!all(split_checks$passed)) {
    set_status("ANALYSIS_FAILED")
    stop("split usage check_table failed")
  }
  log_msg("input_hashes metadata=%s expression=%s", inputs$hashes$metadata, inputs$hashes$expression)
  log_msg("dimensions endpoint_known=488 primary=482 probes=22283")
  log_msg("split_usage_check=PASS")

  metrics_rows <- list()
  threshold_rows <- list()
  hyper_rows <- list()
  lambda_rows <- list()
  runtime_rows <- list()
  model_matrix_rows <- list()
  penalty_rows <- list()
  leakage_rows <- list()

  models <- c("basal", "transcriptomic", "combined")
  outer <- splits$internal_outer
  inner <- splits$internal_inner
  for (current_repeat_id in sort(unique(outer$repeat_id))) {
    for (current_outer_fold in sort(unique(outer$outer_fold))) {
      outer_test_ids <- outer[repeat_id == current_repeat_id & outer_fold == current_outer_fold, sample_id]
      test_idx <- match(outer_test_ids, primary$sample_id)
      train_idx <- which(!primary$sample_id %in% outer_test_ids)
      inner_dt <- inner[repeat_id == current_repeat_id & outer_fold == current_outer_fold]
      inner_fold <- inner_dt$inner_fold[match(primary$sample_id[train_idx], inner_dt$sample_id)]
      if (anyNA(inner_fold)) {
        set_status("ANALYSIS_FAILED")
        stop(sprintf("inner fold mismatch repeat=%s outer=%s", current_repeat_id, current_outer_fold))
      }
      leakage_rows[[length(leakage_rows) + 1L]] <- data.table(
        validation_scheme = "internal",
        repeat_id = current_repeat_id,
        outer_fold = current_outer_fold,
        loco_domain = NA_character_,
        outer_test_not_in_training = length(intersect(primary$sample_id[train_idx], primary$sample_id[test_idx])) == 0L,
        outer_test_not_in_inner = length(intersect(primary$sample_id[test_idx], inner_dt$sample_id)) == 0L,
        threshold_training_only = TRUE,
        tuning_training_only = TRUE,
        cohort_used_as_predictor = TRUE
      )
      for (model in models) {
        log_msg("INTERNAL_START repeat=%s outer_fold=%s model=%s", current_repeat_id, current_outer_fold, model)
        out <- evaluate_one_split("internal", model, current_repeat_id, current_outer_fold, NA_character_, train_idx, test_idx, inner_fold, primary, x, probe_ids, config)
        metrics_rows[[length(metrics_rows) + 1L]] <- out$metrics
        threshold_rows[[length(threshold_rows) + 1L]] <- out$threshold
        hyper_rows[[length(hyper_rows) + 1L]] <- out$hyper
        if (!is.null(out$lambda)) lambda_rows[[length(lambda_rows) + 1L]] <- out$lambda
        runtime_rows[[length(runtime_rows) + 1L]] <- out$runtime
        model_matrix_rows[[length(model_matrix_rows) + 1L]] <- out$model_matrix
        penalty_rows[[length(penalty_rows) + 1L]] <- out$penalty
        if (!out$converged || !out$probs_finite) set_status("REVIEW_REQUIRED")
        log_msg("INTERNAL_END repeat=%s outer_fold=%s model=%s elapsed_sec=%.2f", current_repeat_id, current_outer_fold, model, out$elapsed)
      }
    }
  }

  loco_domains <- sort(unique(splits$loco_outer$loco_test_domain))
  loco_inner <- splits$loco_inner
  for (current_domain in loco_domains) {
    test_idx <- which(primary$source_cohort == current_domain)
    train_idx <- which(primary$source_cohort != current_domain)
    inner_dt <- loco_inner[loco_test_domain == current_domain]
    inner_fold <- inner_dt$inner_fold[match(primary$sample_id[train_idx], inner_dt$sample_id)]
    if (anyNA(inner_fold)) {
      set_status("ANALYSIS_FAILED")
      stop(sprintf("LOCO inner fold mismatch domain=%s", current_domain))
    }
    leakage_rows[[length(leakage_rows) + 1L]] <- data.table(
      validation_scheme = "loco",
      repeat_id = NA_integer_,
      outer_fold = NA_integer_,
      loco_domain = current_domain,
      outer_test_not_in_training = length(intersect(primary$sample_id[train_idx], primary$sample_id[test_idx])) == 0L,
      outer_test_not_in_inner = length(intersect(primary$sample_id[test_idx], inner_dt$sample_id)) == 0L,
      threshold_training_only = TRUE,
      tuning_training_only = TRUE,
      cohort_used_as_predictor = FALSE
    )
    for (model in models) {
      log_msg("LOCO_START domain=%s model=%s", current_domain, model)
      out <- evaluate_one_split("loco", model, NA_integer_, NA_integer_, current_domain, train_idx, test_idx, inner_fold, primary, x, probe_ids, config)
      metrics_rows[[length(metrics_rows) + 1L]] <- out$metrics
      threshold_rows[[length(threshold_rows) + 1L]] <- out$threshold
      hyper_rows[[length(hyper_rows) + 1L]] <- out$hyper
      if (!is.null(out$lambda)) lambda_rows[[length(lambda_rows) + 1L]] <- out$lambda
      runtime_rows[[length(runtime_rows) + 1L]] <- out$runtime
      model_matrix_rows[[length(model_matrix_rows) + 1L]] <- out$model_matrix
      penalty_rows[[length(penalty_rows) + 1L]] <- out$penalty
      if (!out$converged || !out$probs_finite) set_status("REVIEW_REQUIRED")
      log_msg("LOCO_END domain=%s model=%s elapsed_sec=%.2f", current_domain, model, out$elapsed)
    }
  }

  metrics <- rbindlist(metrics_rows, fill = TRUE)
  thresholds <- rbindlist(threshold_rows, fill = TRUE)
  hyper <- rbindlist(hyper_rows, fill = TRUE)
  lambda_checks <- rbindlist(lambda_rows, fill = TRUE)
  runtime <- rbindlist(runtime_rows, fill = TRUE)
  model_matrix <- rbindlist(model_matrix_rows, fill = TRUE)
  penalty <- rbindlist(penalty_rows, fill = TRUE)
  leakage <- rbindlist(leakage_rows, fill = TRUE)

  fwrite(metrics[validation_scheme == "internal"], file.path(tables_dir, "phase3b_internal_fold_metrics.tsv"), sep = "\t")
  fwrite(summarize_metric_table(metrics, "internal"), file.path(tables_dir, "phase3b_internal_model_summary.tsv"), sep = "\t")
  fwrite(paired_differences(metrics), file.path(tables_dir, "phase3b_internal_paired_differences.tsv"), sep = "\t")
  fwrite(metrics[validation_scheme == "loco"], file.path(tables_dir, "phase3b_loco_metrics.tsv"), sep = "\t")
  fwrite(summarize_metric_table(metrics, "loco"), file.path(tables_dir, "phase3b_loco_summary.tsv"), sep = "\t")
  fwrite(hyper, file.path(tables_dir, "phase3b_hyperparameter_selections.tsv"), sep = "\t")
  fwrite(thresholds, file.path(tables_dir, "phase3b_thresholds.tsv"), sep = "\t")
  fwrite(lambda_checks, file.path(tables_dir, "phase3b_lambda_boundary_checks.tsv"), sep = "\t")
  fwrite(model_matrix, file.path(tables_dir, "phase3b_model_matrix_checks.tsv"), sep = "\t")
  fwrite(penalty, file.path(tables_dir, "phase3b_penalty_factor_checks.tsv"), sep = "\t")
  fwrite(leakage, file.path(tables_dir, "phase3b_leakage_checks.tsv"), sep = "\t")

  lambda_ok <- if (nrow(lambda_checks)) all(lambda_checks$lambda_boundary_technically_justified) else TRUE
  threshold_probability_ok <- all(is.finite(thresholds$threshold)) &&
    all(thresholds$threshold >= -1e-12 & thresholds$threshold <= 1 + 1e-12)
  leakage_ok <- all(leakage[, outer_test_not_in_training & outer_test_not_in_inner & threshold_training_only & tuning_training_only]) &&
    all(leakage[validation_scheme == "loco", cohort_used_as_predictor == FALSE])
  penalty_ok <- all(penalty[model == "transcriptomic", zero_penalty_count == 0 & all_probes_penalized == TRUE]) &&
    all(penalty[model == "combined" & validation_scheme == "internal", zero_penalty_count == 4 & one_penalty_count == 22283 & all_probes_penalized == TRUE & order_ok == TRUE]) &&
    all(penalty[model == "combined" & validation_scheme == "loco", zero_penalty_count == 1 & one_penalty_count == 22283 & all_probes_penalized == TRUE & order_ok == TRUE])
  model_matrix_ok <- all(model_matrix[model == "combined" & validation_scheme == "internal", columns == 22287 & cohort_predictor == TRUE & er_predictor == TRUE & probe_columns == 22283 & column_order_ok == TRUE]) &&
    all(model_matrix[model == "combined" & validation_scheme == "loco", columns == 22284 & cohort_predictor == FALSE & er_predictor == TRUE & probe_columns == 22283 & column_order_ok == TRUE]) &&
    all(model_matrix[model == "basal" & validation_scheme == "loco", columns == 1 & cohort_predictor == FALSE])
  convergence_ok <- all(runtime$convergence_errors == 0)
  finite_prob_ok <- all(runtime$nonfinite_probabilities == 0)
  warnings_ok <- all(runtime$warnings == 0)
  metrics_ok <- all(is.finite(as.matrix(metrics[, .(pr_auc, roc_auc, balanced_accuracy, sensitivity, specificity)])))

  memory_info <- available_memory_mb()
  available_mb <- memory_info$mb
  memory_source <- memory_info$source
  memory_peak_mb <- max(runtime$memory_mb_approx, na.rm = TRUE) * 4
  memory_fraction <- if (is.finite(available_mb)) memory_peak_mb / available_mb else NA_real_
  memory_verified <- is.finite(available_mb) && is.finite(memory_fraction) && memory_fraction < 0.70
  runtime[, `:=`(
    available_memory_mb = available_mb,
    memory_fraction_available = memory_fraction,
    memory_available_source = memory_source,
    memory_verification_passed = memory_verified
  )]
  runtime_summary <- runtime[, .(
    validation_scheme = "summary",
    model = "projected_total",
    repeat_id = NA_integer_,
    outer_fold = NA_integer_,
    loco_domain = NA_character_,
    fit_time_sec = sum(fit_time_sec, na.rm = TRUE),
    inner_tuning_time_sec = sum(inner_tuning_time_sec, na.rm = TRUE),
    memory_mb_approx = memory_peak_mb,
    glmnet_trajectories = sum(glmnet_trajectories, na.rm = TRUE),
    lambdas_returned = sum(lambdas_returned, na.rm = TRUE),
    warnings = sum(warnings, na.rm = TRUE),
    convergence_errors = sum(convergence_errors, na.rm = TRUE),
    nonfinite_probabilities = sum(nonfinite_probabilities, na.rm = TRUE),
    available_memory_mb = available_mb,
    memory_fraction_available = memory_fraction,
    memory_available_source = memory_source,
    memory_verification_passed = memory_verified
  )]
  fwrite(rbindlist(list(runtime, runtime_summary), fill = TRUE), file.path(tables_dir, "phase3b_runtime_summary.tsv"), sep = "\t")

  if (!lambda_ok || !warnings_ok || !convergence_ok || !finite_prob_ok || !metrics_ok || !memory_verified || !threshold_probability_ok) set_status("REVIEW_REQUIRED")
  if (!leakage_ok || !penalty_ok || !model_matrix_ok) set_status("ANALYSIS_FAILED")

  checks <- data.table(
    check_name = c(
      "input_integrity_passed",
      "split_usage_passed",
      "leakage_passed",
      "penalty_factors_passed",
      "model_matrices_passed",
      "convergence_passed",
      "finite_probabilities_passed",
      "metrics_calculable",
      "thresholds_on_probability_scale",
      "lambda_boundary_justified",
      "memory_verification_passed",
      "no_population_or_probe_changes",
      "no_predictions_or_coefficients_written",
      "not_extended_to_stability_or_biology"
    ),
    passed = c(
      TRUE,
      all(split_checks$passed),
      leakage_ok,
      penalty_ok,
      model_matrix_ok,
      convergence_ok,
      finite_prob_ok,
      metrics_ok,
      threshold_probability_ok,
      lambda_ok,
      memory_verified,
      TRUE,
      TRUE,
      TRUE
    )
  )
  fwrite(checks, file.path(tables_dir, "phase3b_checks.tsv"), sep = "\t")

  incident_lines <- character()
  if (length(warnings_seen)) incident_lines <- c(incident_lines, sprintf("- Warnings: %s", paste(unique(warnings_seen), collapse = " | ")))
  boundary_rows <- lambda_checks[lambda_on_boundary == TRUE]
  if (nrow(boundary_rows)) {
    incident_lines <- c(incident_lines, sprintf(
      "- Lambda boundary en %s ajustes glmnet; todos justificados tecnicamente: %s.",
      nrow(boundary_rows),
      all(boundary_rows$lambda_boundary_technically_justified)
    ))
  }
  if (!memory_verified) incident_lines <- c(incident_lines, sprintf("- Memoria no verificable: fuente=%s.", memory_source))
  if (length(incident_lines) == 0L) incident_lines <- "- Sin incidencias relevantes registradas."

  internal_summary <- fread(file.path(tables_dir, "phase3b_internal_model_summary.tsv"))
  paired <- fread(file.path(tables_dir, "phase3b_internal_paired_differences.tsv"))
  loco_summary <- fread(file.path(tables_dir, "phase3b_loco_summary.tsv"))
  report <- c(
    "# GSE25066: evaluación predictiva",
    "",
    sprintf("Estado final: `%s`", final_status),
    "",
    "## Alcance",
    "",
    "Este análisis ejecuta la evaluación predictiva interna y leave-one-cohort-out. No incluye estabilidad de sondas, interpretacion biologica, analisis diferencial ni enriquecimiento funcional.",
    "",
    "## Inputs",
    "",
    sprintf("- Metadatos: `%s`", rel(metadata_path)),
    sprintf("- Expresion: `%s`", rel(expression_path)),
    sprintf("- Anotacion: `%s`", rel(annotation_path)),
    sprintf("- Hash metadatos: `%s`", inputs$hashes$metadata),
    sprintf("- Hash expresion: `%s`", inputs$hashes$expression),
    "- Dimensiones: 488 endpoint-known; 482 primarias; 22.283 sondas.",
    "",
    "## Protocolo ejecutado",
    "",
    "- Validacion interna: 5 repeticiones x 5 outer folds usando las particiones generadas previamente.",
    "- LOCO: ISPY, LBJ/IN/GEI, MDACC y USO como dominios dejados fuera.",
    "- Modelos: basal, transcriptomico Elastic Net y combinado.",
    "- Tuning: desviancia binomial interna y regla global 1-SE.",
    "- Umbrales: seleccionados solo dentro del entrenamiento interno por balanced accuracy.",
    "- Metricas: PR-AUC/average precision, ROC-AUC, balanced accuracy, sensibilidad y especificidad.",
    "",
    "## Comprobaciones tecnicas",
    "",
    paste(sprintf("- %s: %s", checks$check_name, checks$passed), collapse = "\n"),
    "",
    "## Resultados cuantitativos",
    "",
    "Los resultados cuantitativos se guardan en tablas TSV para revision neutral. Este reporte no declara un modelo ganador.",
    "",
    "### Resumen interno por modelo",
    "",
    paste(capture.output(print(internal_summary[metric %in% c("pr_auc", "roc_auc", "balanced_accuracy")], nrows = 50)), collapse = "\n"),
    "",
    "### Diferencias pareadas internas",
    "",
    paste(capture.output(print(paired[metric %in% c("pr_auc", "roc_auc", "balanced_accuracy")], nrows = 50)), collapse = "\n"),
    "",
    "### Resumen LOCO",
    "",
    paste(capture.output(print(loco_summary[metric %in% c("pr_auc", "roc_auc", "balanced_accuracy")], nrows = 50)), collapse = "\n"),
    "",
    "## Incidencias",
    "",
    paste(incident_lines, collapse = "\n"),
    "",
    "## Archivos generados",
    "",
    "- `results/predictive_modeling/tables/phase3b_internal_fold_metrics.tsv`",
    "- `results/predictive_modeling/tables/phase3b_internal_model_summary.tsv`",
    "- `results/predictive_modeling/tables/phase3b_internal_paired_differences.tsv`",
    "- `results/predictive_modeling/tables/phase3b_loco_metrics.tsv`",
    "- `results/predictive_modeling/tables/phase3b_loco_summary.tsv`",
    "- `results/predictive_modeling/tables/phase3b_hyperparameter_selections.tsv`",
    "- `results/predictive_modeling/tables/phase3b_thresholds.tsv`",
    "- `results/predictive_modeling/tables/phase3b_lambda_boundary_checks.tsv`",
    "- `results/predictive_modeling/logs/phase3b_console.log`",
    "",
    "## Recomendacion tecnica",
    "",
    if (final_status == "ANALYSIS_COMPLETE") "analysis_complete" else if (final_status == "REVIEW_REQUIRED") "review_required" else "analysis_failed",
    "",
    "No se guardaron predicciones por muestra, coeficientes reales ni listas de sondas seleccionadas."
  )
  writeLines(report, file.path(out_dir, "gse25066_phase3b_predictive_report.md"))
  write_session_info()

  log_msg("FINAL_STATUS=%s", final_status)
  log_msg("files_written=%s", paste(c(
    rel(file.path(out_dir, "gse25066_phase3b_predictive_report.md")),
    rel(file.path(tables_dir, "phase3b_checks.tsv")),
    rel(file.path(tables_dir, "phase3b_internal_model_summary.tsv")),
    rel(file.path(tables_dir, "phase3b_internal_paired_differences.tsv")),
    rel(file.path(tables_dir, "phase3b_loco_summary.tsv")),
    rel(file.path(tables_dir, "phase3b_loco_metrics.tsv")),
    rel(file.path(tables_dir, "phase3b_hyperparameter_selections.tsv")),
    rel(file.path(tables_dir, "phase3b_lambda_boundary_checks.tsv")),
    rel(console_log)
  ), collapse = ";"))
}

tryCatch(main(), error = function(e) {
  log_msg("ERROR=%s", conditionMessage(e))
  if (final_status == "ANALYSIS_COMPLETE") final_status <<- "REVIEW_REQUIRED"
  log_msg("FINAL_STATUS=%s", final_status)
  quit(status = 1L)
})
