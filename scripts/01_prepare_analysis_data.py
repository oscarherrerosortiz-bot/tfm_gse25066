#!/usr/bin/env python3
"""Prepare GSE25066 analysis data for the predictive TFM workflow.

This script performs data preparation only. It does not run EDA, PCA,
feature selection, differential expression, predictive modelling, enrichment, scaling,
batch correction or imputation.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import platform
import re
import sys
import urllib.request
from collections import Counter
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "GSE25066"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "GSE25066"
RESULTS_DIR = PROJECT_ROOT / "results" / "data_preparation"

SERIES_MATRIX = RAW_DIR / "GSE25066_series_matrix.txt.gz"
CLEAN_METADATA = PROCESSED_DIR / "gse25066_sample_metadata_clean.tsv"
GPL96_ANNOT = RAW_DIR / "GPL96.annot.gz"
GPL96_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"

ANALYSIS_METADATA_OUT = PROCESSED_DIR / "gse25066_analysis_metadata.tsv"
EXPRESSION_ENDPOINT_OUT = PROCESSED_DIR / "gse25066_expression_probe_endpoint_known.tsv.gz"
PROBE_ANNOTATION_OUT = PROCESSED_DIR / "gpl96_probe_annotation.tsv.gz"
PREP_REPORT_OUT = RESULTS_DIR / "gse25066_data_preparation_report.md"
EXCLUSIONS_OUT = RESULTS_DIR / "gse25066_sample_exclusions.tsv"

EXPECTED = {
    "n_expression_samples": 508,
    "n_expression_probes": 22283,
    "n_endpoint_known": 488,
    "n_endpoint_pcr": 99,
    "n_endpoint_rd": 389,
    "n_primary": 482,
    "n_primary_pcr": 98,
    "n_primary_rd": 384,
    "n_er_indeterminate_endpoint": 4,
    "n_er_missing_endpoint": 2,
}

ALLOWED_RESPONSE_RAW = {"", "pCR", "RD"}
ALLOWED_RESPONSE = {"pCR", "RD"}
ALLOWED_SOURCE = {"ISPY", "LBJ/IN/GEI", "MDACC", "USO"}
ALLOWED_ER_RAW = {"", "positive", "negative", "indeterminate"}
ALLOWED_ER = {"positive", "negative", "indeterminate", "missing"}
MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "unknown", "---"}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_endpoint_expression_hash(
    sample_ids: list[str],
    probe_ids: list[str],
    expression_rows: dict[str, list[str]],
) -> str:
    digest = hashlib.sha256()
    digest.update(("sample_id\t" + "\t".join(probe_ids) + "\n").encode("utf-8"))
    for sample_id in sample_ids:
        digest.update((sample_id + "\t" + "\t".join(expression_rows[sample_id]) + "\n").encode("utf-8"))
    return digest.hexdigest()


def require_inputs() -> None:
    missing = [path for path in (SERIES_MATRIX, CLEAN_METADATA) if not path.exists()]
    if missing:
        formatted = ", ".join(rel(path) for path in missing)
        raise FileNotFoundError(f"Missing required input(s): {formatted}")


def validate_gzip(path: Path) -> None:
    try:
        with gzip.open(path, "rb") as handle:
            handle.read(1024)
    except OSError as exc:
        raise RuntimeError(f"Cannot open gzip file {rel(path)}: {exc}") from exc


def ensure_gpl96_annotation() -> str:
    if GPL96_ANNOT.exists():
        validate_gzip(GPL96_ANNOT)
        return "existing_file"

    tmp_path = GPL96_ANNOT.with_suffix(GPL96_ANNOT.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    try:
        with urllib.request.urlopen(GPL96_URL, timeout=60) as response, tmp_path.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
        validate_gzip(tmp_path)
        tmp_path.replace(GPL96_ANNOT)
        return "downloaded_from_geo"
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(f"Failed to download or validate official GPL96 annotation from {GPL96_URL}.")


def unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def normalize_er(value: object) -> str:
    if value is None:
        return "missing"
    text = str(value).strip()
    if text == "" or text.lower() in MISSING_TOKENS:
        return "missing"
    if text in {"positive", "negative", "indeterminate"}:
        return text
    raise ValueError(f"Unexpected ER category: {text!r}")


def normalize_response(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text in ALLOWED_RESPONSE:
        return text
    if text == "" or text.lower() in MISSING_TOKENS:
        return ""
    raise ValueError(f"Unexpected response category: {text!r}")


def normalize_source(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if text not in ALLOWED_SOURCE:
        raise ValueError(f"Unexpected source_cohort category: {text!r}")
    return text


def read_clean_metadata() -> tuple[list[dict[str, str]], dict[str, list[str]]]:
    with CLEAN_METADATA.open("r", encoding="utf-8", newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))

    geo_ids = [row["geo_accession"] for row in metadata]
    duplicated_geo_ids = [sample for sample, count in Counter(geo_ids).items() if count > 1]
    if duplicated_geo_ids:
        raise RuntimeError(f"Duplicated metadata geo_accession values: {duplicated_geo_ids[:5]}")

    observed = {
        "response_pcr_vs_rd_raw": sorted({row.get("response_pcr_vs_rd", "") for row in metadata}),
        "source_cohort_raw": sorted({row.get("source_cohort", "") for row in metadata}),
        "er_status_raw": sorted({row.get("er_status", "") for row in metadata}),
    }
    if not set(observed["response_pcr_vs_rd_raw"]).issubset(ALLOWED_RESPONSE_RAW):
        raise RuntimeError(f"Unexpected response categories observed: {observed['response_pcr_vs_rd_raw']}")
    if not set(observed["source_cohort_raw"]).issubset(ALLOWED_SOURCE):
        raise RuntimeError(f"Unexpected source_cohort categories observed: {observed['source_cohort_raw']}")
    if not set(observed["er_status_raw"]).issubset(ALLOWED_ER_RAW):
        raise RuntimeError(f"Unexpected er_status categories observed: {observed['er_status_raw']}")

    normalized_rows: list[dict[str, str]] = []
    for row in metadata:
        normalized = dict(row)
        normalized["geo_accession_original"] = row["geo_accession"]
        normalized["original_sample_id"] = row["sample_id"]
        normalized["sample_id"] = row["geo_accession"]
        normalized["response_pcr_vs_rd"] = normalize_response(row.get("response_pcr_vs_rd", ""))
        normalized["source_cohort"] = normalize_source(row.get("source_cohort", ""))
        normalized["er_status"] = normalize_er(row.get("er_status", ""))
        normalized["response_binary"] = {"pCR": "1", "RD": "0"}.get(
            normalized["response_pcr_vs_rd"], ""
        )
        normalized_rows.append(normalized)
    return normalized_rows, observed


def parse_series_matrix() -> tuple[list[str], list[str], dict[str, object], dict[str, list[str]], dict[str, float]]:
    sample_ids: list[str] = []
    probe_ids: list[str] = []
    expression_rows: dict[str, list[str]] = {}
    missing_cells = 0
    non_finite_cells = 0
    total_cells = 0
    min_value = math.inf
    max_value = -math.inf
    numeric_sample: list[float] = []

    with gzip.open(SERIES_MATRIX, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.rstrip("\n") == "!series_matrix_table_begin":
                break
        else:
            raise RuntimeError("Could not find !series_matrix_table_begin in series matrix.")

        header = next(handle).rstrip("\n").split("\t")
        if unquote(header[0]) != "ID_REF":
            raise RuntimeError(f"Unexpected expression header first column: {header[0]!r}")
        sample_ids = [unquote(value) for value in header[1:]]
        if len(sample_ids) != len(set(sample_ids)):
            raise RuntimeError("Duplicated expression sample identifiers detected.")
        expression_rows = {sample_id: [] for sample_id in sample_ids}

        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_end":
                break
            if not line or line.startswith("!"):
                continue
            parts = line.split("\t")
            probe_id = unquote(parts[0])
            values = [unquote(value) for value in parts[1:]]
            if len(values) != len(sample_ids):
                raise RuntimeError(
                    f"Probe {probe_id} has {len(values)} values; expected {len(sample_ids)}."
                )
            probe_ids.append(probe_id)
            for sample_id, value in zip(sample_ids, values):
                text = value.strip()
                total_cells += 1
                if text == "" or text.lower() in MISSING_TOKENS:
                    missing_cells += 1
                    parsed = math.nan
                else:
                    try:
                        parsed = float(text)
                    except ValueError as exc:
                        raise RuntimeError(f"Non-numeric expression value for {probe_id}/{sample_id}: {text!r}") from exc
                    if not math.isfinite(parsed):
                        non_finite_cells += 1
                    else:
                        min_value = min(min_value, parsed)
                        max_value = max(max_value, parsed)
                        if len(numeric_sample) < 250000:
                            numeric_sample.append(parsed)
                expression_rows[sample_id].append(text)

    if len(probe_ids) != len(set(probe_ids)):
        raise RuntimeError("Duplicated expression probe identifiers detected.")

    metrics: dict[str, object] = {
        "n_expression_samples": len(sample_ids),
        "n_expression_probes": len(probe_ids),
        "expression_total_cells": total_cells,
        "expression_missing_cells": missing_cells,
        "expression_missing_pct": 100 * missing_cells / total_cells if total_cells else math.nan,
        "expression_non_finite_cells": non_finite_cells,
        "expression_min": min_value,
        "expression_max": max_value,
    }
    if numeric_sample:
        sorted_sample = sorted(numeric_sample)

        def quantile(probability: float) -> float:
            if len(sorted_sample) == 1:
                return sorted_sample[0]
            position = probability * (len(sorted_sample) - 1)
            lower = math.floor(position)
            upper = math.ceil(position)
            if lower == upper:
                return sorted_sample[int(position)]
            fraction = position - lower
            return sorted_sample[lower] * (1 - fraction) + sorted_sample[upper] * fraction

        metrics.update(
            {
                "expression_quantile_sample_min": quantile(0.0),
                "expression_quantile_sample_q25": quantile(0.25),
                "expression_quantile_sample_median": quantile(0.5),
                "expression_quantile_sample_q75": quantile(0.75),
                "expression_quantile_sample_max": quantile(1.0),
            }
        )
    return sample_ids, probe_ids, metrics, expression_rows, {k: float(v) for k, v in metrics.items() if isinstance(v, (int, float))}


def build_analysis_metadata(metadata: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    output_rows = []
    exclusion_rows = []

    for row in metadata:
        response = row["response_pcr_vs_rd"]
        er_status = row["er_status"]
        endpoint_known = response in ALLOWED_RESPONSE

        if not endpoint_known:
            exclusion_rows.append(
                {
                    "sample_id": row["sample_id"],
                    "exclusion_stage": "endpoint_known",
                    "exclusion_reason": "endpoint_unknown",
                    "response_pcr_vs_rd": "missing",
                    "er_status": er_status,
                    "source_cohort": row["source_cohort"],
                }
            )
            continue

        primary_eligible = er_status in {"positive", "negative"}
        if primary_eligible:
            primary_reason = ""
        elif er_status == "indeterminate":
            primary_reason = "er_indeterminate"
            exclusion_rows.append(
                {
                    "sample_id": row["sample_id"],
                    "exclusion_stage": "primary_modeling",
                    "exclusion_reason": "er_indeterminate",
                    "response_pcr_vs_rd": response,
                    "er_status": er_status,
                    "source_cohort": row["source_cohort"],
                }
            )
        elif er_status == "missing":
            primary_reason = "er_missing"
            exclusion_rows.append(
                {
                    "sample_id": row["sample_id"],
                    "exclusion_stage": "primary_modeling",
                    "exclusion_reason": "er_missing",
                    "response_pcr_vs_rd": response,
                    "er_status": er_status,
                    "source_cohort": row["source_cohort"],
                }
            )
        else:
            raise RuntimeError(f"Unexpected normalized ER category after endpoint filtering: {er_status!r}")

        record = dict(row)
        record["endpoint_known_eligible"] = "True"
        record["primary_modeling_eligible"] = "True" if primary_eligible else "False"
        record["primary_exclusion_reason"] = primary_reason
        output_rows.append(record)

    return output_rows, exclusion_rows


def write_expression_endpoint_matrix(
    analysis_metadata: list[dict[str, str]],
    probe_ids: list[str],
    expression_rows: dict[str, list[str]],
) -> None:
    endpoint_samples = [row["sample_id"] for row in analysis_metadata]
    with gzip.open(EXPRESSION_ENDPOINT_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id"] + probe_ids)
        for sample_id in endpoint_samples:
            writer.writerow([sample_id] + expression_rows[sample_id])


def validate_written_expression_matrix(
    path: Path,
    expected_sample_ids: list[str],
    expected_probe_ids: list[str],
    expression_rows: dict[str, list[str]],
    expected_canonical_hash: str,
) -> dict[str, object]:
    observed_hash = hashlib.sha256()
    n_rows = 0
    observed_sample_ids: list[str] = []
    first_error = ""

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        observed_hash.update(("\t".join(header) + "\n").encode("utf-8"))
        observed_probe_ids = header[1:]

        for row_index, row in enumerate(reader):
            n_rows += 1
            if len(row) != len(header):
                first_error = first_error or f"row {row_index + 2} has {len(row)} columns; expected {len(header)}"
            sample_id = row[0]
            observed_sample_ids.append(sample_id)
            observed_hash.update(("\t".join(row) + "\n").encode("utf-8"))

            if sample_id in expression_rows:
                expected_values = expression_rows[sample_id]
                if row[1:] != expected_values:
                    first_error = first_error or f"expression values differ for sample {sample_id}"
            for value in row[1:]:
                if value == "":
                    first_error = first_error or f"empty expression cell in sample {sample_id}"
                    continue
                try:
                    parsed = float(value)
                except ValueError:
                    first_error = first_error or f"non-numeric expression value in sample {sample_id}: {value!r}"
                    continue
                if not math.isfinite(parsed):
                    first_error = first_error or f"non-finite expression value in sample {sample_id}: {value!r}"

    duplicate_samples = [sample_id for sample_id, count in Counter(observed_sample_ids).items() if count > 1]
    sample_order_ok = observed_sample_ids == expected_sample_ids
    probe_order_ok = observed_probe_ids == expected_probe_ids
    hash_observed = observed_hash.hexdigest()
    hash_match = hash_observed == expected_canonical_hash

    metrics = {
        "derived_matrix_rows": n_rows,
        "derived_matrix_columns": len(header),
        "derived_matrix_probe_columns": len(observed_probe_ids),
        "derived_matrix_duplicate_sample_ids": len(duplicate_samples),
        "derived_matrix_sample_order_ok": sample_order_ok,
        "derived_matrix_probe_order_ok": probe_order_ok,
        "derived_matrix_canonical_hash_expected": expected_canonical_hash,
        "derived_matrix_canonical_hash_observed": hash_observed,
        "derived_matrix_canonical_hash_match": hash_match,
        "derived_matrix_file_sha256": sha256_file(path),
        "derived_matrix_validation_error": first_error,
    }
    if n_rows != len(expected_sample_ids):
        raise RuntimeError(f"Written expression matrix has {n_rows} rows; expected {len(expected_sample_ids)}.")
    if len(header) != len(expected_probe_ids) + 1:
        raise RuntimeError(f"Written expression matrix has {len(header)} columns; expected {len(expected_probe_ids) + 1}.")
    if duplicate_samples:
        raise RuntimeError(f"Written expression matrix has duplicated sample IDs: {duplicate_samples[:5]}.")
    if not sample_order_ok:
        raise RuntimeError("Written expression matrix sample order differs from analysis metadata.")
    if not probe_order_ok:
        raise RuntimeError("Written expression matrix probe order differs from series matrix.")
    if not hash_match:
        raise RuntimeError("Written expression matrix canonical hash differs from expected hash.")
    if first_error:
        raise RuntimeError(f"Written expression matrix validation failed: {first_error}")
    return metrics


def normalize_annotation_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def get_column_name(columns: list[str], candidates: list[str]) -> str | None:
    normalized = {normalize_annotation_key(column): column for column in columns}
    for candidate in candidates:
        key = normalize_annotation_key(candidate)
        if key in normalized:
            return normalized[key]
    return None


def split_gene_symbols(raw_value: str) -> list[str]:
    text = (raw_value or "").strip()
    if text == "" or text.lower() in MISSING_TOKENS:
        return []
    pieces = re.split(r"\s*///\s*", text)
    symbols: list[str] = []
    for piece in pieces:
        symbol = piece.strip()
        if not symbol or symbol.lower() in MISSING_TOKENS:
            continue
        if symbol not in symbols:
            symbols.append(symbol)
    return symbols


def read_gpl96_annotation(probe_ids: list[str]) -> tuple[list[dict[str, object]], dict[str, object], dict[str, str]]:
    table_header: list[str] | None = None
    annotation_by_probe: dict[str, dict[str, str]] = {}
    duplicate_annotation_ids: list[str] = []
    in_table = False

    with gzip.open(GPL96_ANNOT, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                table_header = next(handle).rstrip("\n").split("\t")
                continue
            if line == "!platform_table_end":
                break
            if not in_table or not table_header or line.startswith("!") or not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(table_header):
                parts += [""] * (len(table_header) - len(parts))
            row = dict(zip(table_header, parts))
            probe_id = row.get("ID", "")
            if probe_id:
                if probe_id in annotation_by_probe:
                    duplicate_annotation_ids.append(probe_id)
                    continue
                annotation_by_probe[probe_id] = row

    if not table_header:
        raise RuntimeError("Could not find GPL96 platform table in official annotation.")
    if duplicate_annotation_ids:
        examples = ", ".join(duplicate_annotation_ids[:10])
        raise RuntimeError(f"Duplicated ID values found in GPL96 annotation table: {examples}")

    symbol_col = get_column_name(table_header, ["Gene Symbol", "Gene symbol", "GeneSymbols"])
    title_col = get_column_name(table_header, ["Gene Title", "Gene title"])
    gene_id_col = get_column_name(table_header, ["Gene ID", "Entrez Gene", "Entrez Gene ID"])
    unigene_col = get_column_name(table_header, ["UniGene ID", "Unigene ID"])
    genbank_col = get_column_name(table_header, ["GenBank Accession", "GenBank Accession Number"])
    refseq_col = get_column_name(table_header, ["RefSeq Transcript ID", "RefSeq Accession"])

    if symbol_col is None:
        raise RuntimeError("GPL96 annotation does not contain a recognizable Gene Symbol column.")

    records = []
    symbols_across_all_mappings = set()
    symbols_from_single_symbol_probes = set()
    unmatched_probe_ids: list[str] = []
    for probe_id in probe_ids:
        annotation_record_found = probe_id in annotation_by_probe
        if not annotation_record_found:
            unmatched_probe_ids.append(probe_id)
        row = annotation_by_probe.get(probe_id, {})
        gene_symbol_raw = row.get(symbol_col, "")
        symbols = split_gene_symbols(gene_symbol_raw)
        symbols_across_all_mappings.update(symbols)
        n_symbols = len(symbols)
        if n_symbols == 0:
            mapping_status = "no_symbol"
            single_symbol = ""
        elif n_symbols == 1:
            mapping_status = "single_symbol"
            single_symbol = symbols[0]
            symbols_from_single_symbol_probes.add(single_symbol)
        else:
            mapping_status = "multiple_symbols"
            single_symbol = ""

        records.append(
            {
                "probe_id": probe_id,
                "annotation_record_found": "True" if annotation_record_found else "False",
                "gene_symbol_raw": gene_symbol_raw,
                "gene_title_raw": row.get(title_col, "") if title_col else "",
                "gene_id_raw": row.get(gene_id_col, "") if gene_id_col else "",
                "unigene_id_raw": row.get(unigene_col, "") if unigene_col else "",
                "genbank_accession_raw": row.get(genbank_col, "") if genbank_col else "",
                "refseq_transcript_id_raw": row.get(refseq_col, "") if refseq_col else "",
                "gene_symbols_normalized": " /// ".join(symbols),
                "n_gene_symbols": n_symbols,
                "mapping_status": mapping_status,
                "gene_symbol_single": single_symbol,
            }
        )

    status_counts = Counter(str(row["mapping_status"]) for row in records)
    with_symbol = sum(int(row["n_gene_symbols"]) > 0 for row in records)
    metrics = {
        "n_gpl96_annotation_ids": len(annotation_by_probe),
        "annotation_rows_total": len(records),
        "annotation_probe_count": len(records),
        "n_expression_probes_matched_to_gpl96": len(probe_ids) - len(unmatched_probe_ids),
        "n_expression_probes_unmatched_to_gpl96": len(unmatched_probe_ids),
        "expression_probes_unmatched_to_gpl96": ", ".join(unmatched_probe_ids),
        "n_duplicate_gpl96_ids": 0,
        "probes_with_symbol": with_symbol,
        "probes_without_symbol": int(status_counts["no_symbol"]),
        "probes_single_symbol": int(status_counts["single_symbol"]),
        "probes_multiple_symbols": int(status_counts["multiple_symbols"]),
        "unique_symbols_across_all_mappings": len(symbols_across_all_mappings),
        "unique_symbols_from_single_symbol_probes": len(symbols_from_single_symbol_probes),
        "symbol_coverage_pct": round(100 * with_symbol / len(records), 2),
    }
    columns_used = {
        "symbol_col": symbol_col or "",
        "title_col": title_col or "",
        "gene_id_col": gene_id_col or "",
        "unigene_col": unigene_col or "",
        "genbank_col": genbank_col or "",
        "refseq_col": refseq_col or "",
    }
    return records, metrics, columns_used


def check_expected_counts(metrics: dict[str, int | float]) -> list[str]:
    issues = []
    for key, expected in EXPECTED.items():
        observed = metrics.get(key)
        if observed != expected:
            issues.append(f"{key}: observed {observed}, expected {expected}")
    return issues


def check_annotation_expectations(annotation_metrics: dict[str, object]) -> list[str]:
    issues = []
    with_symbol = int(annotation_metrics["probes_with_symbol"])
    multi = int(annotation_metrics["probes_multiple_symbols"])
    coverage = float(annotation_metrics["symbol_coverage_pct"])
    matched = int(annotation_metrics["n_expression_probes_matched_to_gpl96"])
    unmatched = int(annotation_metrics["n_expression_probes_unmatched_to_gpl96"])
    all_symbols = int(annotation_metrics["unique_symbols_across_all_mappings"])
    single_symbols = int(annotation_metrics["unique_symbols_from_single_symbol_probes"])
    if matched != 22283:
        issues.append(f"Expression probes matched to GPL96 differs from expected: {matched} vs 22283")
    if unmatched != 0:
        issues.append(
            "Expression probes unmatched to GPL96: "
            + str(annotation_metrics["expression_probes_unmatched_to_gpl96"])
        )
    if abs(with_symbol - 21156) > 25:
        issues.append(f"GPL96 probes with symbol differs from prior check: {with_symbol} vs around 21156")
    if abs(multi - 1223) > 25:
        issues.append(f"GPL96 multiple-symbol probes differs from prior check: {multi} vs around 1223")
    if abs(all_symbols - 13909) > 25:
        issues.append(
            f"Unique symbols across all mappings differs from expected: {all_symbols} vs around 13909"
        )
    if abs(single_symbols - 12502) > 25:
        issues.append(
            "Unique symbols from single-symbol probes differs from expected: "
            f"{single_symbols} vs around 12502"
        )
    if not (94.5 <= coverage <= 95.3):
        issues.append(f"GPL96 symbol coverage outside expected range: {coverage}% vs around 94.94%")
    return issues


def write_report(
    status: str,
    checksums: dict[str, str],
    expression_metrics: dict[str, object],
    counts: dict[str, int],
    response_counts: Counter,
    primary_counts: Counter,
    er_counts: Counter,
    cohort_counts: Counter,
    observed_categories: dict[str, list[str]],
    annotation_metrics: dict[str, object],
    annotation_columns: dict[str, str],
    derived_matrix_metrics: dict[str, object],
    exclusion_counts: Counter,
    issues: list[str],
    outputs: list[Path],
) -> None:
    endpoint_prevalence = 100 * counts["n_endpoint_pcr"] / counts["n_endpoint_known"]
    primary_prevalence = 100 * counts["n_primary_pcr"] / counts["n_primary"]
    lines = [
        "# GSE25066 data preparation report",
        "",
        f"## 1. Estado final",
        "",
        f"**{status}**",
        "",
        "## 2. Inputs y checksums",
        "",
        f"- `{rel(SERIES_MATRIX)}`: `{checksums['series_matrix_sha256']}`",
        f"- `{rel(CLEAN_METADATA)}`: `{checksums['clean_metadata_sha256']}`",
        f"- `{rel(GPL96_ANNOT)}`: `{checksums['gpl96_annotation_sha256']}`",
        "",
        "## 3. Dimensiones",
        "",
        f"- Expresion original: {expression_metrics['n_expression_samples']} muestras x {expression_metrics['n_expression_probes']} sondas.",
        f"- Missingness de expresion: {expression_metrics['expression_missing_cells']} celdas ({expression_metrics['expression_missing_pct']:.6f}%).",
        f"- Valores no finitos: {expression_metrics['expression_non_finite_cells']}.",
        f"- Rango global de expresion: min={expression_metrics['expression_min']}, max={expression_metrics['expression_max']}.",
        f"- Endpoint-known: {counts['n_endpoint_known']} muestras.",
        f"- Conjunto primario: {counts['n_primary']} muestras.",
        f"- Anotacion: {annotation_metrics['annotation_probe_count']} sondas.",
        "",
        "## 4. Poblacion endpoint-known",
        "",
        f"- n = {counts['n_endpoint_known']}",
        f"- pCR = {response_counts.get('pCR', 0)}",
        f"- RD = {response_counts.get('RD', 0)}",
        f"- Prevalencia pCR = {endpoint_prevalence:.2f}%",
        "",
        "## 5. Poblacion primaria de modelado",
        "",
        f"- n = {counts['n_primary']}",
        f"- pCR = {primary_counts.get('pCR', 0)}",
        f"- RD = {primary_counts.get('RD', 0)}",
        f"- Prevalencia pCR = {primary_prevalence:.2f}%",
        "",
        "## 6. Distribucion de ER y cohorte",
        "",
        "- ER: " + ", ".join(f"{key}={value}" for key, value in sorted(er_counts.items())),
        "- Cohorte: " + ", ".join(f"{key}={value}" for key, value in sorted(cohort_counts.items())),
        "",
        "## 7. Correspondencia expresion-metadatos",
        "",
        "- Todas las muestras de expresion tienen metadatos.",
        "- No hay metadatos duplicados por identificador GEO.",
        "- La matriz endpoint-known se ha escrito en el mismo orden que `gse25066_analysis_metadata.tsv`.",
        "",
        "## 8. Integridad de la matriz derivada",
        "",
        f"- Filas de muestras verificadas tras reabrir: {derived_matrix_metrics['derived_matrix_rows']}",
        f"- Columnas verificadas tras reabrir: {derived_matrix_metrics['derived_matrix_columns']}",
        f"- Columnas de sondas: {derived_matrix_metrics['derived_matrix_probe_columns']}",
        f"- IDs de muestra duplicados: {derived_matrix_metrics['derived_matrix_duplicate_sample_ids']}",
        f"- Orden de muestras correcto: {derived_matrix_metrics['derived_matrix_sample_order_ok']}",
        f"- Orden de sondas correcto: {derived_matrix_metrics['derived_matrix_probe_order_ok']}",
        f"- Hash canonico esperado: `{derived_matrix_metrics['derived_matrix_canonical_hash_expected']}`",
        f"- Hash canonico del archivo generado: `{derived_matrix_metrics['derived_matrix_canonical_hash_observed']}`",
        f"- Igualdad de hashes canonicos: {derived_matrix_metrics['derived_matrix_canonical_hash_match']}",
        f"- SHA-256 del archivo comprimido generado: `{derived_matrix_metrics['derived_matrix_file_sha256']}`",
        "- Igualdad exacta del contenido respecto al series matrix: confirmada.",
        "",
        "## 9. Categorias observadas y normalizadas",
        "",
        f"- Respuesta observada en metadatos limpios: {observed_categories['response_pcr_vs_rd_raw']}.",
        "- Respuesta normalizada: `pCR`, `RD`; blancos tratados como endpoint desconocido.",
        f"- Cohorte observada: {observed_categories['source_cohort_raw']}.",
        "- Cohorte normalizada: `ISPY`, `LBJ/IN/GEI`, `MDACC`, `USO`.",
        f"- ER observado: {observed_categories['er_status_raw']}.",
        "- ER normalizado: `positive`, `negative`, `indeterminate`, `missing`.",
        "",
        "## 10. Exclusiones por etapa",
        "",
    ]
    for key, value in sorted(exclusion_counts.items()):
        lines.append(f"- {key}: {value}")
    lines += [
        "",
        "Las muestras con ER indeterminado o missing no se denominan excluidas del estudio; siguen en el conjunto endpoint-known y solo quedan fuera del modelado primario.",
        "",
        "## 11. Integridad de la anotacion",
        "",
        f"- IDs oficiales GPL96: {annotation_metrics['n_gpl96_annotation_ids']}",
        f"- Sondas de expresion encontradas en GPL96: {annotation_metrics['n_expression_probes_matched_to_gpl96']}",
        f"- Sondas de expresion no encontradas en GPL96: {annotation_metrics['n_expression_probes_unmatched_to_gpl96']}",
        f"- Registros duplicados en GPL96: {annotation_metrics['n_duplicate_gpl96_ids']}",
        f"- Sondas con simbolo: {annotation_metrics['probes_with_symbol']}",
        f"- Sondas sin simbolo: {annotation_metrics['probes_without_symbol']}",
        f"- Sondas con simbolo unico: {annotation_metrics['probes_single_symbol']}",
        f"- Sondas con varios simbolos: {annotation_metrics['probes_multiple_symbols']}",
        f"- Simbolos distintos considerando todos los mapeos: {annotation_metrics['unique_symbols_across_all_mappings']}",
        f"- Simbolos distintos procedentes de sondas con mapeo univoco: {annotation_metrics['unique_symbols_from_single_symbol_probes']}",
        f"- Cobertura con simbolo: {annotation_metrics['symbol_coverage_pct']}%",
        f"- Columnas oficiales usadas: {annotation_columns}",
        "",
        "## 12. Archivos generados",
        "",
    ]
    lines += [f"- `{rel(path)}`" for path in outputs]
    lines += [
        "",
        "## 13. Advertencias o discrepancias",
        "",
    ]
    if issues:
        lines += [f"- {issue}" for issue in issues]
    else:
        lines.append("- No se detectaron discrepancias frente a los checks preespecificados.")
    lines += [
        "",
        "## 14. Analisis no realizados",
        "",
        "No se han realizado EDA, PCA, filtrado por varianza, estandarizacion, batch correction, analisis diferencial ni modelado.",
    ]
    PREP_REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plain_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_gzip_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_plain_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_existing_tsv_matches(path: Path, expected_rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    if not path.exists():
        raise RuntimeError(f"Expected existing file is absent and will not be recreated in this QC run: {rel(path)}")
    observed_rows = read_plain_tsv(path)
    normalized_expected = [
        {field: str(row.get(field, "")) for field in fieldnames}
        for row in expected_rows
    ]
    normalized_observed = [
        {field: str(row.get(field, "")) for field in fieldnames}
        for row in observed_rows
    ]
    if normalized_observed != normalized_expected:
        raise RuntimeError(
            f"Existing file {rel(path)} differs from regenerated content; not overwriting."
        )



def main() -> int:
    print("Iniciando la preparación de datos para GSE25066.")
    print(f"Input series matrix: {rel(SERIES_MATRIX)}")
    print(f"Input clean metadata: {rel(CLEAN_METADATA)}")

    require_inputs()
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    series_sha_before = sha256_file(SERIES_MATRIX)
    metadata_sha = sha256_file(CLEAN_METADATA)

    annotation_mode = ensure_gpl96_annotation()
    gpl96_sha = sha256_file(GPL96_ANNOT)

    metadata_rows, observed_categories = read_clean_metadata()
    expression_sample_ids, probe_ids, expression_metrics, expression_rows, _ = parse_series_matrix()

    if expression_metrics["n_expression_samples"] != EXPECTED["n_expression_samples"]:
        raise RuntimeError("Unexpected number of expression samples.")
    if expression_metrics["n_expression_probes"] != EXPECTED["n_expression_probes"]:
        raise RuntimeError("Unexpected number of expression probes.")
    if expression_metrics["expression_missing_cells"] != 0:
        raise RuntimeError("Expression matrix contains missing values.")
    if expression_metrics["expression_non_finite_cells"] != 0:
        raise RuntimeError("Expression matrix contains non-finite values.")

    expression_set = set(expression_sample_ids)
    metadata_set = {row["sample_id"] for row in metadata_rows}
    expression_without_metadata = sorted(expression_set - metadata_set)
    metadata_without_expression = sorted(metadata_set - expression_set)
    if expression_without_metadata:
        raise RuntimeError(f"Expression samples without metadata: {expression_without_metadata[:10]}")
    if metadata_without_expression:
        raise RuntimeError(f"Metadata samples without expression: {metadata_without_expression[:10]}")

    analysis_metadata, exclusions = build_analysis_metadata(metadata_rows)
    endpoint_sample_ids = [row["sample_id"] for row in analysis_metadata]
    if any(sample_id not in expression_rows for sample_id in endpoint_sample_ids):
        raise RuntimeError("At least one endpoint-known sample is missing from expression rows.")

    response_counts = Counter(row["response_pcr_vs_rd"] for row in analysis_metadata)
    primary_metadata = [row for row in analysis_metadata if row["primary_modeling_eligible"] == "True"]
    primary_counts = Counter(row["response_pcr_vs_rd"] for row in primary_metadata)
    er_counts_endpoint = Counter(row["er_status"] for row in analysis_metadata)
    cohort_counts_endpoint = Counter(row["source_cohort"] for row in analysis_metadata)
    exclusion_counts = Counter(row["exclusion_reason"] for row in exclusions)

    counts = {
        "n_endpoint_known": len(analysis_metadata),
        "n_endpoint_pcr": response_counts.get("pCR", 0),
        "n_endpoint_rd": response_counts.get("RD", 0),
        "n_primary": len(primary_metadata),
        "n_primary_pcr": primary_counts.get("pCR", 0),
        "n_primary_rd": primary_counts.get("RD", 0),
        "n_er_indeterminate_endpoint": er_counts_endpoint.get("indeterminate", 0),
        "n_er_missing_endpoint": er_counts_endpoint.get("missing", 0),
    }

    count_issues = check_expected_counts(
        {
            "n_expression_samples": expression_metrics["n_expression_samples"],
            "n_expression_probes": expression_metrics["n_expression_probes"],
            **counts,
        }
    )

    endpoint_sample_ids = [row["sample_id"] for row in analysis_metadata]
    expected_matrix_hash = canonical_endpoint_expression_hash(endpoint_sample_ids, probe_ids, expression_rows)
    write_expression_endpoint_matrix(analysis_metadata, probe_ids, expression_rows)
    derived_matrix_metrics = validate_written_expression_matrix(
        EXPRESSION_ENDPOINT_OUT,
        endpoint_sample_ids,
        probe_ids,
        expression_rows,
        expected_matrix_hash,
    )

    analysis_fields = list(analysis_metadata[0].keys())
    exclusion_fields = [
        "sample_id",
        "exclusion_stage",
        "exclusion_reason",
        "response_pcr_vs_rd",
        "er_status",
        "source_cohort",
    ]
    write_plain_tsv(ANALYSIS_METADATA_OUT, analysis_metadata, analysis_fields)
    write_plain_tsv(EXCLUSIONS_OUT, exclusions, exclusion_fields)
    verify_existing_tsv_matches(ANALYSIS_METADATA_OUT, analysis_metadata, analysis_fields)
    verify_existing_tsv_matches(EXCLUSIONS_OUT, exclusions, exclusion_fields)

    annotation_rows, annotation_metrics, annotation_columns = read_gpl96_annotation(probe_ids)
    write_gzip_tsv(
        PROBE_ANNOTATION_OUT,
        annotation_rows,
        [
            "probe_id",
            "annotation_record_found",
            "gene_symbol_raw",
            "gene_title_raw",
            "gene_id_raw",
            "unigene_id_raw",
            "genbank_accession_raw",
            "refseq_transcript_id_raw",
            "gene_symbols_normalized",
            "n_gene_symbols",
            "mapping_status",
            "gene_symbol_single",
        ],
    )
    annotation_issues = check_annotation_expectations(annotation_metrics)

    series_sha_after = sha256_file(SERIES_MATRIX)
    if series_sha_after != series_sha_before:
        raise RuntimeError("Raw series matrix checksum changed during execution.")

    checksums = {
        "series_matrix_sha256": series_sha_after,
        "clean_metadata_sha256": metadata_sha,
        "gpl96_annotation_sha256": gpl96_sha,
    }

    issues = count_issues + annotation_issues
    status = "REVIEW_REQUIRED" if issues else "ANALYSIS_COMPLETE"
    outputs = [
        ANALYSIS_METADATA_OUT,
        EXCLUSIONS_OUT,
        EXPRESSION_ENDPOINT_OUT,
        PROBE_ANNOTATION_OUT,
        PREP_REPORT_OUT,
    ]

    write_report(
        status,
        checksums,
        expression_metrics,
        counts,
        response_counts,
        primary_counts,
        er_counts_endpoint,
        cohort_counts_endpoint,
        observed_categories,
        annotation_metrics,
        annotation_columns,
        derived_matrix_metrics,
        exclusion_counts,
        issues,
        outputs,
    )
    print(f"Annotation mode: {annotation_mode}")
    print(
        f"Original dimensions: {expression_metrics['n_expression_samples']} samples x "
        f"{expression_metrics['n_expression_probes']} probes"
    )
    print("Expression-metadata correspondence: OK")
    print(
        f"Endpoint-known: n={counts['n_endpoint_known']}, "
        f"pCR={counts['n_endpoint_pcr']}, RD={counts['n_endpoint_rd']}"
    )
    print(
        f"Primary modeling: n={counts['n_primary']}, "
        f"pCR={counts['n_primary_pcr']}, RD={counts['n_primary_rd']}"
    )
    print(
        "Exclusions: "
        + ", ".join(f"{key}={value}" for key, value in sorted(exclusion_counts.items()))
    )
    print(
        "Derived matrix verified: "
        f"rows={derived_matrix_metrics['derived_matrix_rows']}, "
        f"columns={derived_matrix_metrics['derived_matrix_columns']}, "
        f"probe_columns={derived_matrix_metrics['derived_matrix_probe_columns']}"
    )
    print(f"Sample order matches metadata: {derived_matrix_metrics['derived_matrix_sample_order_ok']}")
    print(f"Probe order matches series matrix: {derived_matrix_metrics['derived_matrix_probe_order_ok']}")
    print(
        "Canonical expression hash comparison: "
        f"expected={derived_matrix_metrics['derived_matrix_canonical_hash_expected']}, "
        f"observed={derived_matrix_metrics['derived_matrix_canonical_hash_observed']}, "
        f"match={derived_matrix_metrics['derived_matrix_canonical_hash_match']}"
    )
    print(
        "GPL96 records: "
        f"official_ids={annotation_metrics['n_gpl96_annotation_ids']}, "
        f"matched={annotation_metrics['n_expression_probes_matched_to_gpl96']}, "
        f"unmatched={annotation_metrics['n_expression_probes_unmatched_to_gpl96']}, "
        f"duplicate_ids={annotation_metrics['n_duplicate_gpl96_ids']}"
    )
    print(
        "Probe annotation mapping: "
        f"single_symbol={annotation_metrics['probes_single_symbol']}, "
        f"multiple_symbols={annotation_metrics['probes_multiple_symbols']}, "
        f"no_symbol={annotation_metrics['probes_without_symbol']}"
    )
    print(
        "Distinct symbols: "
        f"all_mappings={annotation_metrics['unique_symbols_across_all_mappings']}, "
        f"single_symbol_probes={annotation_metrics['unique_symbols_from_single_symbol_probes']}"
    )
    print("Outputs created:")
    for path in outputs:
        print(f"- {rel(path)}")
    print(f"Final state: {status}")
    print(f"STATUS: {status}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
