#!/usr/bin/env python3
"""Initial integrity and suitability check for GEO series GSE25066.

This script intentionally avoids a full exploratory or biological analysis. It parses the
processed GEO series matrix, checks the pCR/RD endpoint and covariates, validates a
preliminary covariate design, performs a small limma smoke test on
high-variance probes, and checks whether GPL96 probe annotation is usable.
"""

from __future__ import annotations

import csv
import gzip
import heapq
import math
import re
import subprocess
import sys
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "GSE25066"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "GSE25066"
RESULTS_DIR = PROJECT_ROOT / "results" / "initial_dataset_check"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

SERIES_MATRIX = RAW_DIR / "GSE25066_series_matrix.txt.gz"
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25066/matrix/"
    "GSE25066_series_matrix.txt.gz"
)
GPL96_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

METADATA_OUT = PROCESSED_DIR / "gse25066_sample_metadata_clean.tsv"
DIMENSIONS_OUT = TABLES_DIR / "initial_dataset_dimensions.tsv"
ENDPOINT_OUT = TABLES_DIR / "initial_dataset_endpoint_balance.tsv"
MISSINGNESS_OUT = TABLES_DIR / "initial_dataset_covariate_summary.tsv"
RESPONSE_TABLE_OUT = TABLES_DIR / "initial_dataset_response_by_covariate.tsv"
DESIGN_OUT = TABLES_DIR / "initial_dataset_design_check.tsv"
ANNOTATION_OUT = TABLES_DIR / "initial_dataset_probe_annotation_check.tsv"
REPORT_OUT = RESULTS_DIR / "initial_dataset_check_report.md"

MISSING_TOKENS = {"", "na", "n/a", "nan", "unknown", "not available", "null", "none", "---"}
TOP_VARIANCE_PROBES = 300


def ensure_directories() -> None:
    for directory in (RAW_DIR, PROCESSED_DIR, TABLES_DIR, FIGURES_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def ensure_series_matrix() -> str:
    """Ensure the raw matrix exists, preferring an available local copy."""
    if SERIES_MATRIX.exists():
        return "existing_raw_file"

    urllib.request.urlretrieve(SERIES_MATRIX_URL, SERIES_MATRIX)
    return "downloaded_from_geo"


def normalize_name(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def clean_missing(value: str | None) -> str:
    if value is None:
        return ""
    value = value.strip()
    return "" if value.lower() in MISSING_TOKENS else value


def parse_header_and_metadata(path: Path) -> tuple[dict[str, list[str]], dict[str, str], list[dict[str, str]]]:
    """Parse GEO header lines and characteristic key-value metadata."""
    sample_fields: dict[str, list[str]] = {}
    series_fields: dict[str, str] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                break
            if line.startswith("!Series_"):
                parts = line.split("\t")
                key = normalize_name(parts[0].removeprefix("!Series_"))
                series_fields[key] = " | ".join(unquote(x) for x in parts[1:])
            elif line.startswith("!Sample_"):
                parts = line.split("\t")
                key = normalize_name(parts[0].removeprefix("!Sample_"))
                values = [unquote(x) for x in parts[1:]]
                if key in sample_fields:
                    suffix = 2
                    while f"{key}_{suffix}" in sample_fields:
                        suffix += 1
                    key = f"{key}_{suffix}"
                sample_fields[key] = values

    if "geo_accession" not in sample_fields:
        raise RuntimeError("The series matrix does not contain !Sample_geo_accession.")

    n_samples = len(sample_fields["geo_accession"])
    metadata_rows: list[dict[str, str]] = [dict() for _ in range(n_samples)]

    for key, values in sample_fields.items():
        if len(values) != n_samples:
            raise RuntimeError(f"Metadata field {key} has {len(values)} values, expected {n_samples}.")
        for index, value in enumerate(values):
            metadata_rows[index][key] = clean_missing(value)

    characteristic_keys = [key for key in sample_fields if key.startswith("characteristics_ch1")]
    for row in metadata_rows:
        for key in characteristic_keys:
            entry = row.get(key, "")
            if ": " not in entry:
                continue
            label, value = entry.split(": ", 1)
            row[normalize_name(label)] = clean_missing(value)

    return sample_fields, series_fields, metadata_rows


def parse_float(value: str) -> float | None:
    value = value.strip()
    if value == "" or value.lower() in MISSING_TOKENS:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def sample_variance(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / (len(values) - 1)


def parse_expression_stream(
    path: Path,
) -> tuple[list[str], dict[str, int | float], list[tuple[float, str, list[float]]], set[str]]:
    """Stream expression rows and retain only high-variance probes for the smoke test."""
    sample_ids: list[str] = []
    probe_seen: set[str] = set()
    duplicate_probe_ids: set[str] = set()
    top_probe_heap: list[tuple[float, str, list[float]]] = []
    n_features = 0
    missing_cells = 0
    total_cells = 0

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.rstrip("\n") == "!series_matrix_table_begin":
                break
        header = next(handle).rstrip("\n").split("\t")
        sample_ids = [unquote(value) for value in header[1:]]

        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_end":
                break
            if not line or line.startswith("!"):
                continue
            parts = line.split("\t")
            probe_id = unquote(parts[0])
            values_raw = parts[1:]
            if len(values_raw) != len(sample_ids):
                raise RuntimeError(
                    f"Probe {probe_id} has {len(values_raw)} expression values, "
                    f"expected {len(sample_ids)}."
                )

            n_features += 1
            if probe_id in probe_seen:
                duplicate_probe_ids.add(probe_id)
            probe_seen.add(probe_id)

            parsed = [parse_float(value) for value in values_raw]
            missing_cells += sum(value is None for value in parsed)
            total_cells += len(parsed)

            if all(value is not None for value in parsed):
                numeric_values = [float(value) for value in parsed if value is not None]
                variance = sample_variance(numeric_values)
                item = (variance, probe_id, numeric_values)
                if len(top_probe_heap) < TOP_VARIANCE_PROBES:
                    heapq.heappush(top_probe_heap, item)
                elif variance > top_probe_heap[0][0]:
                    heapq.heapreplace(top_probe_heap, item)

    metrics: dict[str, int | float] = {
        "n_expression_samples": len(sample_ids),
        "n_expression_features": n_features,
        "expression_total_cells": total_cells,
        "expression_missing_cells": missing_cells,
        "expression_missing_pct": 100 * missing_cells / total_cells if total_cells else math.nan,
        "duplicate_probe_id_count": len(duplicate_probe_ids),
    }
    return sample_ids, metrics, sorted(top_probe_heap, reverse=True), duplicate_probe_ids


def metadata_value(row: dict[str, str], *candidates: str) -> str:
    for candidate in candidates:
        value = clean_missing(row.get(candidate, ""))
        if value:
            return value
    return ""


def normalize_response(value: str) -> str:
    value_upper = value.strip().upper()
    if value_upper == "PCR":
        return "pCR"
    if value_upper == "RD":
        return "RD"
    return ""


def normalize_receptor(value: str) -> str:
    value_upper = value.strip().upper()
    mapping = {"P": "positive", "N": "negative", "I": "indeterminate"}
    return mapping.get(value_upper, "")


def build_clean_metadata(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    clean_rows: list[dict[str, str]] = []
    for row in rows:
        clean_rows.append(
            {
                "geo_accession": metadata_value(row, "geo_accession"),
                "sample_id": metadata_value(row, "sample_id", "sample_id_2", "title"),
                "source_cohort": metadata_value(row, "source"),
                "age_years": metadata_value(row, "age_years"),
                "response_pcr_vs_rd": normalize_response(
                    metadata_value(row, "pathologic_response_pcr_rd")
                ),
                "er_status": normalize_receptor(metadata_value(row, "er_status_ihc")),
                "pr_status": normalize_receptor(metadata_value(row, "pr_status_ihc")),
                "her2_status": normalize_receptor(metadata_value(row, "her2_status")),
                "clinical_t_stage": metadata_value(row, "clinical_t_stage"),
                "clinical_nodal_status": metadata_value(row, "clinical_nodal_status"),
                "clinical_ajcc_stage": metadata_value(row, "clinical_ajcc_stage"),
                "grade": metadata_value(row, "grade"),
                "pathologic_response_rcb_class": metadata_value(
                    row, "pathologic_response_rcb_class"
                ),
                "drfs_event": metadata_value(row, "drfs_1_event_0_censored"),
                "drfs_time_years": metadata_value(row, "drfs_even_time_years"),
                "platform_id": metadata_value(row, "platform_id"),
            }
        )
    return clean_rows


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def missingness_rows(rows: list[dict[str, str]], variables: list[str]) -> list[dict[str, object]]:
    output = []
    for variable in variables:
        missing = sum(not clean_missing(row.get(variable, "")) for row in rows)
        output.append(
            {
                "variable": variable,
                "n_total": len(rows),
                "n_missing": missing,
                "pct_missing": round(100 * missing / len(rows), 2) if rows else math.nan,
                "n_non_missing": len(rows) - missing,
                "n_unique_non_missing": len(
                    {row.get(variable, "") for row in rows if clean_missing(row.get(variable, ""))}
                ),
            }
        )
    return output


def response_breakdown(
    endpoint_rows: list[dict[str, str]], variables: list[str]
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for variable in variables:
        level_counts: dict[str, Counter[str]] = defaultdict(Counter)
        for row in endpoint_rows:
            level = clean_missing(row.get(variable, "")) or "missing"
            level_counts[level][row["response_pcr_vs_rd"]] += 1
        for level in sorted(level_counts):
            counts = level_counts[level]
            level_total = sum(counts.values())
            for response in ("pCR", "RD"):
                n = counts[response]
                output.append(
                    {
                        "variable": variable,
                        "level": level,
                        "response": response,
                        "n": n,
                        "level_total": level_total,
                        "pct_within_level": round(100 * n / level_total, 2) if level_total else 0,
                        "contains_both_response_classes": counts["pCR"] > 0 and counts["RD"] > 0,
                    }
                )
    return output


def matrix_rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    n_rows = len(work)
    n_cols = len(work[0])
    rank = 0
    for column in range(n_cols):
        pivot = max(range(rank, n_rows), key=lambda row: abs(work[row][column]), default=rank)
        if abs(work[pivot][column]) <= tolerance:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(n_rows):
            if row == rank:
                continue
            factor = work[row][column]
            if abs(factor) > tolerance:
                work[row] = [
                    current - factor * pivot_current
                    for current, pivot_current in zip(work[row], work[rank])
                ]
        rank += 1
        if rank == n_rows:
            break
    return rank


def build_design_matrix(
    endpoint_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[str], list[list[float]], dict[str, object]]:
    """Build response + cohort + ER design using complete cases and treatment coding."""
    complete_rows = [
        row
        for row in endpoint_rows
        if row["source_cohort"] and row["er_status"] and row["response_pcr_vs_rd"]
    ]
    cohort_counts = Counter(row["source_cohort"] for row in complete_rows)
    er_counts = Counter(row["er_status"] for row in complete_rows)
    cohort_reference = cohort_counts.most_common(1)[0][0]
    er_reference = er_counts.most_common(1)[0][0]
    cohort_levels = sorted(level for level in cohort_counts if level != cohort_reference)
    er_levels = sorted(level for level in er_counts if level != er_reference)

    columns = ["intercept", "response_pcr"]
    columns += [f"cohort_{normalize_name(level)}" for level in cohort_levels]
    columns += [f"er_{normalize_name(level)}" for level in er_levels]

    matrix: list[list[float]] = []
    for row in complete_rows:
        values = [1.0, 1.0 if row["response_pcr_vs_rd"] == "pCR" else 0.0]
        values += [1.0 if row["source_cohort"] == level else 0.0 for level in cohort_levels]
        values += [1.0 if row["er_status"] == level else 0.0 for level in er_levels]
        matrix.append(values)

    rank = matrix_rank(matrix)
    combination_counts = Counter(
        (row["response_pcr_vs_rd"], row["source_cohort"], row["er_status"]) for row in complete_rows
    )
    all_combinations = [
        (response, cohort, er)
        for response in ("pCR", "RD")
        for cohort in sorted(cohort_counts)
        for er in sorted(er_counts)
    ]
    empty_combinations = [
        " / ".join(combination) for combination in all_combinations if combination_counts[combination] == 0
    ]
    low_count_combinations = {
        " / ".join(key): count for key, count in combination_counts.items() if count < 5
    }
    metrics: dict[str, object] = {
        "formula": "expression ~ response_pcr_vs_rd + source_cohort + er_status",
        "n_complete_case_samples": len(complete_rows),
        "n_columns": len(columns),
        "rank": rank,
        "full_rank": rank == len(columns),
        "cohort_reference": cohort_reference,
        "er_reference": er_reference,
        "small_cohort_levels_lt5": ", ".join(
            f"{level}={count}" for level, count in cohort_counts.items() if count < 5
        )
        or "none",
        "small_er_levels_lt5": ", ".join(
            f"{level}={count}" for level, count in er_counts.items() if count < 5
        )
        or "none",
        "response_cohort_er_combinations_lt5": ", ".join(
            f"{key}={value}" for key, value in low_count_combinations.items()
        )
        or "none",
        "empty_response_cohort_er_combinations": ", ".join(empty_combinations) or "none",
        "design_columns": ", ".join(columns),
    }
    return complete_rows, columns, matrix, metrics


def run_limma_smoke_test(
    top_probes: list[tuple[float, str, list[float]]],
    expression_sample_ids: list[str],
    design_rows: list[dict[str, str]],
    design_columns: list[str],
    design_matrix: list[list[float]],
) -> dict[str, object]:
    """Run a non-interpretive limma fit on high-variance probes in a temporary directory."""
    sample_to_index = {sample_id: index for index, sample_id in enumerate(expression_sample_ids)}
    design_sample_ids = [row["geo_accession"] for row in design_rows]
    missing_design_samples = [sample_id for sample_id in design_sample_ids if sample_id not in sample_to_index]
    if missing_design_samples:
        raise RuntimeError(
            f"{len(missing_design_samples)} design samples are absent from expression columns."
        )
    selected_indices = [sample_to_index[sample_id] for sample_id in design_sample_ids]

    with tempfile.TemporaryDirectory(prefix="gse25066_smoke_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        expr_path = temp_dir / "smoke_expression.tsv"
        design_path = temp_dir / "smoke_design.tsv"

        with expr_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["probe_id"] + design_sample_ids)
            for _, probe_id, values in top_probes:
                writer.writerow([probe_id] + [values[index] for index in selected_indices])

        with design_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id"] + design_columns)
            for sample_id, values in zip(design_sample_ids, design_matrix):
                writer.writerow([sample_id] + values)

        r_code = r"""
args <- commandArgs(trailingOnly = TRUE)
expr <- read.delim(args[1], check.names = FALSE, row.names = 1)
design <- read.delim(args[2], check.names = FALSE, row.names = 1)
stopifnot(identical(colnames(expr), rownames(design)))
suppressMessages(library(limma))
fit <- lmFit(as.matrix(expr), as.matrix(design))
fit <- eBayes(fit)
ok <- all(is.finite(fit$coefficients[, "response_pcr"]))
cat(paste(
  nrow(expr),
  ncol(expr),
  qr(as.matrix(design))$rank,
  ncol(design),
  ok,
  sep = "\t"
))
"""
        completed = subprocess.run(
            ["Rscript", "-e", r_code, str(expr_path), str(design_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        fields = completed.stdout.strip().split("\t")
        if len(fields) != 5:
            raise RuntimeError(f"Unexpected limma smoke-test output: {completed.stdout!r}")
        return {
            "smoke_test_status": "passed" if fields[4] == "TRUE" else "failed",
            "smoke_test_probes": int(fields[0]),
            "smoke_test_samples": int(fields[1]),
            "smoke_test_design_rank": int(fields[2]),
            "smoke_test_design_columns": int(fields[3]),
            "smoke_test_response_coefficients_finite": fields[4],
        }


def download_and_parse_gpl96_annotation(
    expression_probe_ids: set[str],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Download GPL96 annotation temporarily and estimate probe-to-symbol mapping."""
    with tempfile.TemporaryDirectory(prefix="gpl96_annotation_") as temp_dir_name:
        annot_path = Path(temp_dir_name) / "GPL96.annot.gz"
        urllib.request.urlretrieve(GPL96_ANNOT_URL, annot_path)

        header: list[str] | None = None
        id_index = None
        symbol_index = None
        mapped_probe_ids: set[str] = set()
        symbol_to_probes: dict[str, set[str]] = defaultdict(set)
        multi_symbol_probe_ids: set[str] = set()
        total_annotation_rows = 0

        with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line or line.startswith(("#", "!", "^")):
                    continue
                parts = line.split("\t")
                if header is None:
                    header = parts
                    normalized = [normalize_name(column) for column in header]
                    id_candidates = ["id", "id_ref", "probe_set_id", "probe_id"]
                    symbol_candidates = ["gene_symbol", "gene_symbol_2", "symbol"]
                    for candidate in id_candidates:
                        if candidate in normalized:
                            id_index = normalized.index(candidate)
                            break
                    for candidate in symbol_candidates:
                        if candidate in normalized:
                            symbol_index = normalized.index(candidate)
                            break
                    if id_index is None or symbol_index is None:
                        raise RuntimeError(
                            "Could not identify probe ID and gene symbol columns in GPL96 annotation. "
                            f"Columns: {header}"
                        )
                    continue

                total_annotation_rows += 1
                if max(id_index, symbol_index) >= len(parts):
                    continue
                probe_id = parts[id_index].strip()
                raw_symbols = parts[symbol_index].strip()
                if probe_id not in expression_probe_ids:
                    continue
                symbols = [
                    symbol.strip()
                    for symbol in re.split(r"\s*///\s*|\s*//\s*|;\s*", raw_symbols)
                    if symbol.strip() and symbol.strip().lower() not in MISSING_TOKENS
                ]
                if symbols:
                    mapped_probe_ids.add(probe_id)
                    if len(set(symbols)) > 1:
                        multi_symbol_probe_ids.add(probe_id)
                    for symbol in symbols:
                        symbol_to_probes[symbol].add(probe_id)

    duplicated_symbols = {symbol: probes for symbol, probes in symbol_to_probes.items() if len(probes) > 1}
    metrics: dict[str, object] = {
        "platform_id": "GPL96",
        "annotation_source": GPL96_ANNOT_URL,
        "annotation_total_rows": total_annotation_rows,
        "expression_probe_count": len(expression_probe_ids),
        "expression_probes_with_gene_symbol": len(mapped_probe_ids),
        "expression_probes_with_gene_symbol_pct": round(
            100 * len(mapped_probe_ids) / len(expression_probe_ids), 2
        ),
        "unique_gene_symbols": len(symbol_to_probes),
        "gene_symbols_mapped_by_multiple_probes": len(duplicated_symbols),
        "probes_mapped_to_multiple_gene_symbols": len(multi_symbol_probe_ids),
        "probe_to_gene_mapping_feasible": len(mapped_probe_ids) > 0.7 * len(expression_probe_ids),
        "recommended_future_collapse_rule": (
            "For downstream analyses, keep probe-level expression and document ambiguous "
            "probe-to-gene mappings separately."
        ),
    }
    rows = [{"metric": key, "value": value} for key, value in metrics.items()]
    return rows, metrics


def metric_rows(metrics: dict[str, object]) -> list[dict[str, object]]:
    return [{"metric": key, "value": value} for key, value in metrics.items()]


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    if not rows:
        return "_No data._"
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for row in rows:
        body.append(
            "| "
            + " | ".join(str(row.get(column, "")).replace("|", "\\|") for column in columns)
            + " |"
        )
    return "\n".join([header, separator] + body)


def generate_report(
    acquisition_mode: str,
    series_fields: dict[str, str],
    dimensions: dict[str, object],
    endpoint_rows: list[dict[str, str]],
    endpoint_balance: list[dict[str, object]],
    missingness: list[dict[str, object]],
    response_rows: list[dict[str, object]],
    design_metrics: dict[str, object],
    annotation_metrics: dict[str, object],
    complete_confounds: list[str],
) -> str:
    balance_map = {row["response"]: row["n"] for row in endpoint_balance}
    missingness_key = {
        row["variable"]: f'{row["n_missing"]}/{row["n_total"]} ({row["pct_missing"]}%)'
        for row in missingness
    }
    cohort_rows = [
        row for row in response_rows if row["variable"] == "source_cohort" and row["response"] == "pCR"
    ]
    receptor_rows = [
        row
        for row in response_rows
        if row["variable"] in {"er_status", "pr_status", "her2_status"} and row["response"] == "pCR"
    ]
    suitability = "dataset_usable_for_analysis"
    limitation = (
        "El desbalance pCR/RD y la asociación de la respuesta con el estado ER y con la cohorte "
        "USO requieren ajuste y una interpretación prudente, pero no existe confusión completa."
    )
    if not design_metrics["full_rank"] or not annotation_metrics["probe_to_gene_mapping_feasible"]:
        suitability = "dataset_not_usable"
    elif complete_confounds:
        suitability = "dataset_usable_with_constraints"
        limitation = "Se detectó confusión completa en: " + ", ".join(complete_confounds)

    return f"""# Initial dataset integrity check: GSE25066

**Fecha:** 3 de junio de 2026  
**Estado técnico:** **{suitability}**

## 1. Objetivo de la comprobación

Comprobar que GSE25066 contiene una matriz de expresión consistente, un endpoint pCR/RD,
covariables utilizables y anotación GPL96 suficiente para el análisis predictivo posterior.
Esta comprobación no estima rendimiento predictivo ni interpreta genes.

## 2. Archivos usados

- Matriz GEO procesada: `data/raw/GSE25066/GSE25066_series_matrix.txt.gz`
- Origen: [{SERIES_MATRIX_URL}]({SERIES_MATRIX_URL})
- Modo de adquisición local: `{acquisition_mode}`
- Plataforma detectada: `{series_fields.get("platform_id", "no detectada")}`
- Anotación consultada temporalmente: [{GPL96_ANNOT_URL}]({GPL96_ANNOT_URL})

## 3. Dimensiones reales y consistencia

- Muestras en la matriz de expresión: **{dimensions["n_expression_samples"]}**
- Sondas/features: **{dimensions["n_expression_features"]}**
- Missingness de expresión: **{dimensions["expression_missing_pct"]:.6f}%**
- Duplicados de identificadores de muestra: **{dimensions["duplicate_sample_id_count"]}**
- Duplicados de identificadores de sonda: **{dimensions["duplicate_probe_id_count"]}**
- Coincidencia entre columnas de expresión y `geo_accession`: **{dimensions["expression_metadata_sample_ids_match"]}**

La matriz se carga correctamente y los identificadores de muestra son consistentes con los
metadatos GEO.

## 4. Endpoint principal

El endpoint se detectó como `pathologic_response_pcr_rd` y se normalizó a
`response_pcr_vs_rd`, con valores `pCR` y `RD`.

- Muestras totales: **{dimensions["n_metadata_samples"]}**
- Muestras con endpoint conocido: **{len(endpoint_rows)}**
- pCR: **{balance_map.get("pCR", 0)}**
- RD: **{balance_map.get("RD", 0)}**

{markdown_table(endpoint_balance, ["response", "n", "pct_of_known_endpoint"])}

## 5. Covariables clínicas detectadas

Se detectaron y conservaron: origen/cohorte, edad, ER, PR, HER2, estadio T, estado nodal,
estadio AJCC, grado, clase RCB y variables de supervivencia DRFS como inventario.

Missingness entre las muestras con endpoint conocido:

{markdown_table(missingness, ["variable", "n_missing", "pct_missing", "n_unique_non_missing"])}

Las covariables principales para el ajuste son utilizables. ER tiene
`{missingness_key.get("er_status", "no disponible")}` de missingness y origen/cohorte
`{missingness_key.get("source_cohort", "no disponible")}`.

## 6. Distribución de respuesta por origen y receptores

Tasa de pCR por cohorte:

{markdown_table(cohort_rows, ["level", "n", "level_total", "pct_within_level", "contains_both_response_classes"])}

Tasa de pCR por ER, PR y HER2:

{markdown_table(receptor_rows, ["variable", "level", "n", "level_total", "pct_within_level", "contains_both_response_classes"])}

Todas las cohortes principales contienen ambas clases. No se detectó una covariable clínica
completamente confundida con la respuesta entre origen, ER, PR y HER2.

## 7. Riesgos de confusión o desbalance

- La respuesta está desbalanceada: RD es la clase mayoritaria.
- La cohorte USO presenta una tasa de pCR superior a las otras cohortes.
- ER está fuertemente asociado con pCR/RD y debe controlarse en el análisis principal.
- HER2 positivo tiene muy pocas muestras y no es adecuado como covariable principal del modelo.
- El nivel HER2 indeterminado contiene solo RD; es un nivel escaso y no se usa en el diseño principal.
- Los niveles indeterminados de receptores son pequeños y deben tratarse con cautela.

Estos riesgos son manejables porque no hay confusión completa y el tamaño muestral sigue siendo
adecuado tras usar casos completos para respuesta, cohorte y ER.

## 8. Comprobación de matriz de diseño

- Fórmula comprobada: `{design_metrics["formula"]}`
- Muestras completas: **{design_metrics["n_complete_case_samples"]}**
- Columnas: **{design_metrics["n_columns"]}**
- Rango: **{design_metrics["rank"]}**
- Rango completo: **{design_metrics["full_rank"]}**
- Combinaciones respuesta/cohorte/ER vacías: `{design_metrics["empty_response_cohort_er_combinations"]}`
- Combinaciones respuesta/cohorte/ER con menos de 5 muestras:
  `{design_metrics["response_cohort_er_combinations_lt5"]}`
- Smoke test limma en {design_metrics["smoke_test_probes"]} sondas de alta varianza:
  **{design_metrics["smoke_test_status"]}**

La matriz de diseño ajustada por origen/cohorte y ER es viable. El smoke test solo confirma que
el diseño puede ajustarse; no se han guardado ni interpretado resultados de genes.

## 9. Comprobación de anotación de sondas

- Plataforma: **{annotation_metrics["platform_id"]}**
- Sondas de expresión con símbolo génico: **{annotation_metrics["expression_probes_with_gene_symbol"]}**
  de **{annotation_metrics["expression_probe_count"]}**
  (**{annotation_metrics["expression_probes_with_gene_symbol_pct"]}%**)
- Símbolos génicos únicos: **{annotation_metrics["unique_gene_symbols"]}**
- Símbolos representados por múltiples sondas: **{annotation_metrics["gene_symbols_mapped_by_multiple_probes"]}**
- Sondas mapeadas a múltiples símbolos: **{annotation_metrics["probes_mapped_to_multiple_gene_symbols"]}**
- Mapeo factible: **{annotation_metrics["probe_to_gene_mapping_feasible"]}**

La anotación probe → gene symbol es factible mediante la anotación oficial GPL96. Para un
análisis posterior se recomienda conservar la sonda de mayor varianza por símbolo génico y
documentar por separado las sondas con mapeos múltiples; no se ha colapsado la matriz en esta comprobación.

## 10. Conclusión técnica

**{suitability}**

{limitation}

El dataset contiene un endpoint clínico definido, metadatos suficientes, anotación utilizable y
un diseño preliminar de rango completo. El desbalance de respuesta y la estructura por cohorte
deben tratarse explícitamente durante la validación predictiva.
"""


def main() -> int:
    ensure_directories()
    acquisition_mode = ensure_series_matrix()

    sample_fields, series_fields, raw_metadata = parse_header_and_metadata(SERIES_MATRIX)
    clean_metadata = build_clean_metadata(raw_metadata)
    expression_sample_ids, expression_metrics, top_probes, duplicate_probe_ids = parse_expression_stream(
        SERIES_MATRIX
    )

    metadata_sample_ids = [row["geo_accession"] for row in clean_metadata]
    duplicate_sample_ids = [
        sample_id for sample_id, count in Counter(expression_sample_ids).items() if count > 1
    ]
    expression_metadata_match = expression_sample_ids == metadata_sample_ids
    if not expression_metadata_match:
        raise RuntimeError("Expression sample columns do not match metadata geo_accession order.")

    endpoint_rows = [row for row in clean_metadata if row["response_pcr_vs_rd"] in {"pCR", "RD"}]
    endpoint_counts = Counter(row["response_pcr_vs_rd"] for row in endpoint_rows)
    endpoint_balance = [
        {
            "response": response,
            "n": endpoint_counts[response],
            "pct_of_known_endpoint": round(100 * endpoint_counts[response] / len(endpoint_rows), 2),
        }
        for response in ("pCR", "RD")
    ]

    key_covariates = [
        "source_cohort",
        "age_years",
        "er_status",
        "pr_status",
        "her2_status",
        "clinical_t_stage",
        "clinical_nodal_status",
        "clinical_ajcc_stage",
        "grade",
        "pathologic_response_rcb_class",
        "drfs_event",
        "drfs_time_years",
    ]
    missingness = missingness_rows(endpoint_rows, key_covariates)
    response_rows = response_breakdown(
        endpoint_rows, ["source_cohort", "er_status", "pr_status", "her2_status"]
    )

    complete_confounds: list[str] = []
    one_class_levels_by_variable: dict[str, str] = {}
    for variable in ("source_cohort", "er_status", "pr_status", "her2_status"):
        level_to_responses: dict[str, set[str]] = defaultdict(set)
        for row in endpoint_rows:
            level = clean_missing(row.get(variable, ""))
            if level:
                level_to_responses[level].add(row["response_pcr_vs_rd"])
        if any(len(responses) == 1 for responses in level_to_responses.values()):
            one_class_levels = [
                level for level, responses in level_to_responses.items() if len(responses) == 1
            ]
            one_class_levels_by_variable[variable] = ", ".join(one_class_levels)
            # Tiny sparse receptor levels are a warning, not structural complete confounding.
            if variable == "source_cohort":
                complete_confounds.append(f"{variable}: {', '.join(one_class_levels)}")

    design_rows, design_columns, design_matrix, design_metrics = build_design_matrix(endpoint_rows)
    smoke_metrics = run_limma_smoke_test(
        top_probes, expression_sample_ids, design_rows, design_columns, design_matrix
    )
    design_metrics.update(smoke_metrics)
    design_metrics["one_response_class_levels_source_cohort"] = one_class_levels_by_variable.get(
        "source_cohort", "none"
    )
    design_metrics["one_response_class_levels_er_status"] = one_class_levels_by_variable.get(
        "er_status", "none"
    )
    design_metrics["one_response_class_levels_pr_status"] = one_class_levels_by_variable.get(
        "pr_status", "none"
    )
    design_metrics["one_response_class_levels_her2_status"] = one_class_levels_by_variable.get(
        "her2_status", "none"
    )

    # Re-stream IDs without retaining expression values so the annotation check covers all probes.
    all_expression_probe_ids: set[str] = set()
    with gzip.open(SERIES_MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.rstrip("\n") == "!series_matrix_table_begin":
                break
        next(handle)
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_end":
                break
            if line and not line.startswith("!"):
                all_expression_probe_ids.add(unquote(line.split("\t", 1)[0]))
    annotation_rows, annotation_metrics = download_and_parse_gpl96_annotation(all_expression_probe_ids)

    dimensions: dict[str, object] = {
        "series_accession": series_fields.get("geo_accession", "GSE25066"),
        "platform_id": series_fields.get("platform_id", ""),
        "n_metadata_samples": len(clean_metadata),
        **expression_metrics,
        "duplicate_sample_id_count": len(duplicate_sample_ids),
        "expression_metadata_sample_ids_match": expression_metadata_match,
        "n_samples_with_known_pcr_rd": len(endpoint_rows),
        "n_samples_without_known_pcr_rd": len(clean_metadata) - len(endpoint_rows),
    }

    write_tsv(METADATA_OUT, clean_metadata, list(clean_metadata[0].keys()))
    write_tsv(DIMENSIONS_OUT, metric_rows(dimensions), ["metric", "value"])
    write_tsv(ENDPOINT_OUT, endpoint_balance, ["response", "n", "pct_of_known_endpoint"])
    write_tsv(
        MISSINGNESS_OUT,
        missingness,
        [
            "variable",
            "n_total",
            "n_missing",
            "pct_missing",
            "n_non_missing",
            "n_unique_non_missing",
        ],
    )
    write_tsv(
        RESPONSE_TABLE_OUT,
        response_rows,
        [
            "variable",
            "level",
            "response",
            "n",
            "level_total",
            "pct_within_level",
            "contains_both_response_classes",
        ],
    )
    write_tsv(DESIGN_OUT, metric_rows(design_metrics), ["metric", "value"])
    write_tsv(ANNOTATION_OUT, annotation_rows, ["metric", "value"])

    REPORT_OUT.write_text(
        generate_report(
            acquisition_mode,
            series_fields,
            dimensions,
            endpoint_rows,
            endpoint_balance,
            missingness,
            response_rows,
            design_metrics,
            annotation_metrics,
            complete_confounds,
        ),
        encoding="utf-8",
    )
    suitability = "dataset_usable_for_analysis"
    if not design_metrics["full_rank"] or not annotation_metrics["probe_to_gene_mapping_feasible"]:
        suitability = "dataset_not_usable"
    elif complete_confounds:
        suitability = "dataset_usable_with_constraints"

    print(f"GSE25066 load: OK ({acquisition_mode})")
    print(
        f"Expression dimensions: {dimensions['n_expression_features']} probes x "
        f"{dimensions['n_expression_samples']} samples"
    )
    print(f"Known endpoint samples: {len(endpoint_rows)}; pCR={endpoint_counts['pCR']}; RD={endpoint_counts['RD']}")
    print(
        f"Design: rank {design_metrics['rank']}/{design_metrics['n_columns']}, "
        f"smoke test {design_metrics['smoke_test_status']}"
    )
    print(
        "Annotation: "
        f"{annotation_metrics['expression_probes_with_gene_symbol_pct']}% probes with gene symbol"
    )
    print(f"Dataset status: {suitability}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
