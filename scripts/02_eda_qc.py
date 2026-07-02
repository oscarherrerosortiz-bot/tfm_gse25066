#!/usr/bin/env python3
"""EDA y control de calidad para GSE25066.

This script performs descriptive EDA and sample-level QC before predictive modelling.
It does not train models, perform cross-validation, select probes, run differential
expression, perform enrichment, exclude samples, or alter input data.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import platform
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/gse25066_mpl_config")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/gse25066_xdg_cache")

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import requests
    import scipy
    import sklearn
    from sklearn.decomposition import PCA
except ModuleNotFoundError as exc:
    missing = exc.name or "unknown"
    message = (
        f"Missing required dependency: {missing}. "
        "Create/update the environment from environment.yml and rerun with the same interpreter, "
        "for example: conda env create -f environment.yml && conda activate tfm-gse25066."
    )
    print(message, file=sys.stderr)
    raise SystemExit(2) from exc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
METADATA_PATH = PROJECT_ROOT / "data" / "processed" / "GSE25066" / "gse25066_analysis_metadata.tsv"
EXPRESSION_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "GSE25066"
    / "gse25066_expression_probe_endpoint_known.tsv.gz"
)

RESULTS_DIR = PROJECT_ROOT / "results" / "eda_qc"
TABLES_DIR = RESULTS_DIR
FIGURES_DIR = PROJECT_ROOT / "figures"
COHORT_SUMMARY_OUT = TABLES_DIR / "gse25066_cohort_summary.tsv"
SAMPLE_QC_OUT = TABLES_DIR / "gse25066_sample_qc.tsv"
REPORT_OUT = RESULTS_DIR / "gse25066_eda_qc_report.md"
CLINICAL_FIGURE_OUT = FIGURES_DIR / "gse25066_clinical_overview.png"
PCA_FIGURE_OUT = FIGURES_DIR / "gse25066_pca_structure.png"

RANDOM_STATE = 20260608
PREVIOUS_QC_CANDIDATES = ["GSM615690", "GSM615712", "GSM615746"]
CLINICAL_VARIABLES = [
    "response_pcr_vs_rd",
    "source_cohort",
    "er_status",
    "pr_status",
    "her2_status",
    "age_years",
    "clinical_t_stage",
    "clinical_nodal_status",
    "clinical_ajcc_stage",
    "grade",
]
MISSING_TOKENS = {"", "missing", "na", "n/a", "nan", "none", "null"}
FRIENDLY_VARIABLE_NAMES = {
    "response_pcr_vs_rd": "Respuesta pCR/RD",
    "source_cohort": "Cohorte",
    "er_status": "ER",
    "pr_status": "PR",
    "her2_status": "HER2",
    "age_years": "Edad",
    "clinical_t_stage": "Estadio T",
    "clinical_nodal_status": "Estado nodal",
    "clinical_ajcc_stage": "Estadio AJCC",
    "grade": "Grado",
}
EXPECTED_MISSINGNESS = {
    "er_status": {"missing": 2, "indeterminate": 4},
    "pr_status": {"missing": 2, "indeterminate": 5},
    "her2_status": {"missing": 13, "indeterminate": 4},
    "grade": {"missing": 20, "indeterminate": 15},
}


def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input is missing: {rel(path)}")


def is_missing_value(value: object) -> bool:
    return str(value).strip().lower() in MISSING_TOKENS


def normalized_counter(series: pd.Series) -> Counter:
    normalized = series.astype(str).map(lambda value: "missing" if is_missing_value(value) else value.strip())
    return Counter(normalized)



def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, str]]:
    require_file(METADATA_PATH)
    require_file(EXPRESSION_PATH)

    metadata = pd.read_csv(METADATA_PATH, sep="\t", dtype=str, keep_default_na=False)
    expression = pd.read_csv(EXPRESSION_PATH, sep="\t", compression="gzip", dtype=str, keep_default_na=False)

    if metadata.shape[0] != 488:
        raise RuntimeError(f"Metadata endpoint-known rows differ from expected 488: {metadata.shape[0]}")
    primary_n = int((metadata["primary_modeling_eligible"] == "True").sum())
    if primary_n != 482:
        raise RuntimeError(f"Primary modelling population differs from expected 482: {primary_n}")
    response_counts = Counter(metadata["response_pcr_vs_rd"])
    if response_counts["pCR"] != 99 or response_counts["RD"] != 389:
        raise RuntimeError(f"Endpoint response counts differ from expected 99/389: {response_counts}")
    primary_counts = Counter(metadata.loc[metadata["primary_modeling_eligible"] == "True", "response_pcr_vs_rd"])
    if primary_counts["pCR"] != 98 or primary_counts["RD"] != 384:
        raise RuntimeError(f"Primary response counts differ from expected 98/384: {primary_counts}")

    if expression.shape != (488, 22284):
        raise RuntimeError(f"Expression matrix shape differs from expected 488 x 22284: {expression.shape}")
    if expression["sample_id"].duplicated().any():
        raise RuntimeError("Expression matrix contains duplicated sample IDs.")
    if metadata["sample_id"].duplicated().any():
        raise RuntimeError("Metadata contains duplicated sample IDs.")
    if expression["sample_id"].tolist() != metadata["sample_id"].tolist():
        raise RuntimeError("Expression sample order does not match analysis metadata.")

    probe_ids = expression.columns.tolist()[1:]
    if len(probe_ids) != 22283 or len(set(probe_ids)) != 22283:
        raise RuntimeError("Probe IDs are missing or duplicated in expression matrix.")

    profile_hashes: dict[str, str] = {}
    with gzip.open(EXPRESSION_PATH, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[1:] != probe_ids:
            raise RuntimeError("Probe order changed between pandas read and raw gzip read.")
        for row in reader:
            sample_id = row[0]
            values = row[1:]
            profile_hashes[sample_id] = hashlib.sha256(("\t".join(values)).encode("utf-8")).hexdigest()
            if any(value == "" for value in values):
                raise RuntimeError(f"Expression matrix contains empty cells for sample {sample_id}.")

    expression_numeric = expression[probe_ids].astype(float)
    if not np.isfinite(expression_numeric.to_numpy()).all():
        raise RuntimeError("Expression matrix contains non-finite numeric values.")

    expression_numeric.insert(0, "sample_id", expression["sample_id"].values)
    return metadata, expression_numeric, probe_ids, profile_hashes


def numeric_age(metadata: pd.DataFrame) -> tuple[pd.Series, dict[str, object]]:
    raw_age = metadata["age_years"].map(lambda value: np.nan if is_missing_value(value) else value)
    age = pd.to_numeric(raw_age, errors="coerce")
    non_numeric = int(raw_age.notna().sum() - age.notna().sum())
    finite_mask = np.isfinite(age.to_numpy(dtype=float, na_value=np.nan))
    outside = metadata.loc[age.notna() & ((age < 18) | (age > 100)), "sample_id"].tolist()
    q1 = float(age.quantile(0.25))
    q3 = float(age.quantile(0.75))
    summary = {
        "missing": int(age.isna().sum()),
        "pct_missing": 100 * float(age.isna().mean()),
        "non_numeric": non_numeric,
        "non_finite": int(age.notna().sum() - finite_mask.sum()),
        "min": float(age.min()),
        "max": float(age.max()),
        "median": float(age.median()),
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "outside_18_100_ids": outside,
    }
    return age, summary


def missingness_table(metadata: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for variable in CLINICAL_VARIABLES:
        values = metadata[variable].astype(str)
        missing_mask = values.map(is_missing_value)
        non_missing = values.loc[~missing_mask].map(lambda value: value.strip())
        indeterminate = int(non_missing.str.lower().eq("indeterminate").sum())
        if variable == "grade":
            indeterminate += int(non_missing.eq("4=Indeterminate").sum())
        rows.append(
            {
                "variable": variable,
                "display_name": FRIENDLY_VARIABLE_NAMES.get(variable, variable),
                "n_missing": int(missing_mask.sum()),
                "pct_missing": round(100 * float(missing_mask.mean()), 2),
                "levels": ", ".join(sorted(non_missing.unique().tolist())),
                "n_indeterminate": indeterminate,
            }
        )
    return rows


def q1_q3(series: pd.Series) -> tuple[float, float]:
    return float(series.quantile(0.25)), float(series.quantile(0.75))


def summarize_population(metadata: pd.DataFrame, population: str, source_cohort: str) -> dict[str, object]:
    subset = metadata if source_cohort == "ALL" else metadata.loc[metadata["source_cohort"] == source_cohort]
    age = pd.to_numeric(subset["age_years"].map(lambda value: np.nan if is_missing_value(value) else value), errors="coerce")
    q1, q3 = q1_q3(age)
    response_counts = Counter(subset["response_pcr_vs_rd"])

    def receptor_counts(column: str) -> dict[str, int]:
        counts = normalized_counter(subset[column])
        prefix = column.split("_")[0]
        return {
            f"{prefix}_positive_n": counts["positive"],
            f"{prefix}_negative_n": counts["negative"],
            f"{prefix}_indeterminate_n": counts["indeterminate"],
            f"{prefix}_missing_n": counts["missing"],
        }

    row = {
        "population": population,
        "source_cohort": source_cohort,
        "n": len(subset),
        "pcr_n": response_counts["pCR"],
        "pcr_pct": round(100 * response_counts["pCR"] / len(subset), 2) if len(subset) else 0.0,
        "rd_n": response_counts["RD"],
        "age_median": round(float(age.median()), 3),
        "age_q1": round(q1, 3),
        "age_q3": round(q3, 3),
    }
    row.update(receptor_counts("er_status"))
    row.update(receptor_counts("pr_status"))
    row.update(receptor_counts("her2_status"))
    return row


def cohort_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    endpoint = metadata.copy()
    primary = metadata.loc[metadata["primary_modeling_eligible"] == "True"].copy()
    rows = [summarize_population(endpoint, "endpoint_known", "ALL")]
    rows += [
        summarize_population(endpoint.loc[endpoint["source_cohort"] == cohort], "endpoint_known", cohort)
        for cohort in sorted(endpoint["source_cohort"].unique())
    ]
    rows.append(summarize_population(primary, "primary_modeling", "ALL"))
    rows += [
        summarize_population(primary.loc[primary["source_cohort"] == cohort], "primary_modeling", cohort)
        for cohort in sorted(primary["source_cohort"].unique())
    ]
    return pd.DataFrame(rows)


def robust_z(values: pd.Series, label: str, incidents: list[str]) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    median = numeric.median()
    mad = (numeric - median).abs().median()
    if not np.isfinite(mad) or mad == 0:
        incidents.append(f"MAD zero or non-finite for {label}; robust z-score left as not calculable.")
        return pd.Series([np.nan] * len(numeric), index=values.index)
    return 0.6745 * (numeric - median) / mad


def run_standard_pca(x: np.ndarray) -> tuple[np.ndarray, list[float]]:
    pca = PCA(n_components=10, svd_solver="randomized", random_state=RANDOM_STATE)
    coords = pca.fit_transform(x)
    explained = (100 * pca.explained_variance_ratio_).tolist()
    if coords.shape != (482, 10) or len(explained) != 10:
        raise RuntimeError("PCA did not return expected dimensions.")
    if not np.isfinite(coords).all() or not np.isfinite(explained).all():
        raise RuntimeError("PCA returned non-finite values.")
    return coords, explained


def add_expression_correlations(sample_qc: pd.DataFrame, x: np.ndarray, incidents: list[str]) -> pd.DataFrame:
    corr = np.corrcoef(x)
    if corr.shape != (len(sample_qc), len(sample_qc)):
        raise RuntimeError("Correlation matrix has unexpected dimensions.")
    np.fill_diagonal(corr, np.nan)
    sample_qc["median_correlation_all_samples"] = np.nanmedian(corr, axis=1)

    within_values: list[float] = []
    cohorts = sample_qc["source_cohort"].tolist()
    for idx, cohort in enumerate(cohorts):
        cohort_mask = np.array([other == cohort for other in cohorts])
        cohort_mask[idx] = False
        if not cohort_mask.any():
            within_values.append(np.nan)
        else:
            within_values.append(float(np.nanmedian(corr[idx, cohort_mask])))
    sample_qc["median_correlation_within_cohort"] = within_values
    sample_qc["robust_z_median_correlation_within_cohort"] = np.nan
    for cohort in sorted(sample_qc["source_cohort"].unique()):
        mask = sample_qc["source_cohort"] == cohort
        sample_qc.loc[mask, "robust_z_median_correlation_within_cohort"] = robust_z(
            sample_qc.loc[mask, "median_correlation_within_cohort"],
            f"median_correlation_within_cohort_{cohort}",
            incidents,
        )
    return sample_qc


def run_qc_and_pca(
    metadata: pd.DataFrame,
    expression: pd.DataFrame,
    probe_ids: list[str],
    profile_hashes: dict[str, str],
) -> tuple[pd.DataFrame, list[float], list[str]]:
    incidents: list[str] = []
    primary_metadata = metadata.loc[metadata["primary_modeling_eligible"] == "True"].copy()
    primary_expression = expression.loc[
        expression["sample_id"].isin(primary_metadata["sample_id"]), ["sample_id", *probe_ids]
    ].copy()
    primary_expression = primary_expression.set_index("sample_id").loc[primary_metadata["sample_id"]]
    x = primary_expression.to_numpy(dtype=float)

    coords, explained = run_standard_pca(x)

    sample_qc = primary_metadata[
        ["sample_id", "source_cohort", "response_pcr_vs_rd", "er_status", "age_years"]
    ].copy()
    sample_qc["expression_mean"] = x.mean(axis=1)
    sample_qc["expression_median"] = np.median(x, axis=1)
    sample_qc["expression_sd"] = x.std(axis=1, ddof=1)
    sample_qc["expression_iqr"] = np.percentile(x, 75, axis=1) - np.percentile(x, 25, axis=1)
    sample_qc["robust_z_expression_median"] = robust_z(sample_qc["expression_median"], "expression_median", incidents)
    sample_qc["robust_z_expression_iqr"] = robust_z(sample_qc["expression_iqr"], "expression_iqr", incidents)
    for index in range(5):
        sample_qc[f"PC{index + 1}"] = coords[:, index]
        sample_qc[f"robust_z_PC{index + 1}"] = robust_z(sample_qc[f"PC{index + 1}"], f"PC{index + 1}", incidents)

    sample_qc["expression_profile_hash"] = sample_qc["sample_id"].map(profile_hashes)
    duplicate_hashes = sample_qc["expression_profile_hash"].duplicated(keep=False)
    sample_qc = add_expression_correlations(sample_qc, x, incidents)
    sample_qc["distribution_flag"] = (
        sample_qc["robust_z_expression_median"].abs().gt(5)
        | sample_qc["robust_z_expression_iqr"].abs().gt(5)
    )
    pca_z_cols = [f"robust_z_PC{index}" for index in range(1, 6)]
    sample_qc["pca_extreme_flag"] = sample_qc[pca_z_cols].abs().gt(5).any(axis=1)
    sample_qc["exact_duplicate_profile_flag"] = duplicate_hashes
    sample_qc["qc_candidate_flag"] = (
        sample_qc["distribution_flag"]
        | sample_qc["pca_extreme_flag"]
        | sample_qc["exact_duplicate_profile_flag"]
    )

    numeric_columns = [
        "expression_mean",
        "expression_median",
        "expression_sd",
        "expression_iqr",
        "robust_z_expression_median",
        "robust_z_expression_iqr",
        "PC1",
        "PC2",
        "PC3",
        "PC4",
        "PC5",
        "robust_z_PC1",
        "robust_z_PC2",
        "robust_z_PC3",
        "robust_z_PC4",
        "robust_z_PC5",
        "median_correlation_all_samples",
        "median_correlation_within_cohort",
        "robust_z_median_correlation_within_cohort",
    ]
    for column in numeric_columns:
        sample_qc[column] = sample_qc[column].map(lambda value: "" if pd.isna(value) else round(float(value), 6))
    for column in ["distribution_flag", "pca_extreme_flag", "exact_duplicate_profile_flag", "qc_candidate_flag"]:
        sample_qc[column] = sample_qc[column].map(lambda value: "True" if bool(value) else "False")

    ordered_columns = [
        "sample_id",
        "source_cohort",
        "response_pcr_vs_rd",
        "er_status",
        "age_years",
        "expression_mean",
        "expression_median",
        "expression_sd",
        "expression_iqr",
        "robust_z_expression_median",
        "robust_z_expression_iqr",
        "PC1",
        "PC2",
        "PC3",
        "PC4",
        "PC5",
        "robust_z_PC1",
        "robust_z_PC2",
        "robust_z_PC3",
        "robust_z_PC4",
        "robust_z_PC5",
        "expression_profile_hash",
        "median_correlation_all_samples",
        "median_correlation_within_cohort",
        "robust_z_median_correlation_within_cohort",
        "distribution_flag",
        "pca_extreme_flag",
        "exact_duplicate_profile_flag",
        "qc_candidate_flag",
    ]
    return sample_qc[ordered_columns], explained, incidents


def write_clinical_figure(metadata: pd.DataFrame, age: pd.Series, missing_rows: list[dict[str, object]]) -> None:
    primary = metadata.loc[metadata["primary_modeling_eligible"] == "True"].copy()
    primary_age = age.loc[primary.index]
    rng = np.random.default_rng(RANDOM_STATE)

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 1.05])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])
    colors = {"RD": "#4C78A8", "pCR": "#F58518"}

    cohort_counts = (
        metadata.groupby(["source_cohort", "response_pcr_vs_rd"]).size().unstack(fill_value=0).reindex(columns=["RD", "pCR"])
    )
    bottoms = np.zeros(len(cohort_counts))
    x_pos = np.arange(len(cohort_counts))
    for response in ["RD", "pCR"]:
        ax_a.bar(x_pos, cohort_counts[response].to_numpy(), bottom=bottoms, label=response, color=colors[response])
        for idx, count in enumerate(cohort_counts[response].to_numpy()):
            total = cohort_counts.iloc[idx].sum()
            if count:
                ax_a.text(idx, bottoms[idx] + count / 2, f"{count}\n{100 * count / total:.1f}%", ha="center", va="center", fontsize=8)
        bottoms += cohort_counts[response].to_numpy()
    ax_a.set_title("A. Respuesta por cohorte (n=488)")
    ax_a.set_ylabel("Numero de muestras")
    ax_a.set_xticks(x_pos, cohort_counts.index, rotation=25, ha="right")
    ax_a.legend(title="Respuesta")

    er_subset = primary.loc[primary["er_status"].isin(["negative", "positive"])]
    er_counts = (
        er_subset.groupby(["er_status", "response_pcr_vs_rd"]).size().unstack(fill_value=0).reindex(index=["negative", "positive"], columns=["RD", "pCR"])
    )
    bottoms = np.zeros(len(er_counts))
    x_pos = np.arange(len(er_counts))
    for response in ["RD", "pCR"]:
        ax_b.bar(x_pos, er_counts[response].to_numpy(), bottom=bottoms, label=response, color=colors[response])
        for idx, count in enumerate(er_counts[response].to_numpy()):
            total = er_counts.iloc[idx].sum()
            if count:
                ax_b.text(idx, bottoms[idx] + count / 2, f"{count}\n{100 * count / total:.1f}%", ha="center", va="center", fontsize=8)
        bottoms += er_counts[response].to_numpy()
    ax_b.set_title("B. Respuesta por ER (primaria, n=482)")
    ax_b.set_ylabel("Numero de muestras")
    ax_b.set_xticks(x_pos, ["ER negativo", "ER positivo"])
    ax_b.legend(title="Respuesta")

    age_groups = [
        primary_age.loc[primary["response_pcr_vs_rd"] == response].dropna().astype(float).to_numpy()
        for response in ["RD", "pCR"]
    ]
    ax_c.boxplot(age_groups, tick_labels=["RD", "pCR"], patch_artist=True, boxprops={"facecolor": "#D8DEE9"})
    for idx, values in enumerate(age_groups, start=1):
        jitter = rng.normal(0, 0.04, len(values))
        ax_c.scatter(np.full(len(values), idx) + jitter, values, s=12, alpha=0.55, color="#2F4B7C", linewidths=0)
    ax_c.set_title("C. Edad por respuesta (primaria)")
    ax_c.set_ylabel("Edad (anos)")
    ax_c.set_xlabel("Respuesta")

    missing_df = pd.DataFrame(missing_rows)
    display_order = CLINICAL_VARIABLES
    plot_df = missing_df.set_index("variable").loc[display_order].reset_index()
    y_pos = np.arange(len(plot_df))
    ax_d.barh(y_pos, plot_df["pct_missing"], color="#72B7B2")
    ax_d.set_yticks(y_pos, plot_df["display_name"])
    ax_d.invert_yaxis()
    ax_d.set_xlabel("% missing (denominador: 488 muestras)")
    ax_d.set_title("D. Missingness clinico esencial")
    for idx, row in plot_df.iterrows():
        ax_d.text(row["pct_missing"] + 0.15, idx, f"{int(row['n_missing'])} ({row['pct_missing']:.2f}%)", va="center", fontsize=8)
    ax_d.set_xlim(0, max(6, float(plot_df["pct_missing"].max()) + 2))

    fig.savefig(CLINICAL_FIGURE_OUT, dpi=300)
    plt.close(fig)


def write_pca_figure(sample_qc: pd.DataFrame, explained: list[float]) -> None:
    plot_df = sample_qc.copy()
    plot_df["PC1"] = pd.to_numeric(plot_df["PC1"])
    plot_df["PC2"] = pd.to_numeric(plot_df["PC2"])
    candidates = plot_df.loc[plot_df["qc_candidate_flag"] == "True"]
    panels = [
        ("source_cohort", "A. Cohorte"),
        ("er_status", "B. ER"),
        ("response_pcr_vs_rd", "C. Respuesta"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=True, sharey=True, constrained_layout=True)
    for ax, (column, title) in zip(axes, panels):
        for level in sorted(plot_df[column].unique()):
            subset = plot_df.loc[plot_df[column] == level]
            ax.scatter(subset["PC1"], subset["PC2"], s=18, alpha=0.75, label=level)
        for _, row in candidates.iterrows():
            ax.annotate(row["sample_id"], (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
        ax.set_title(title)
        ax.set_xlabel(f"PC1 ({explained[0]:.2f}% var.)")
        ax.legend(title=column.replace("_", " "), fontsize=8, title_fontsize=8)
    axes[0].set_ylabel(f"PC2 ({explained[1]:.2f}% var.)")
    fig.suptitle("PCA descriptiva sobre 482 muestras primarias")
    fig.savefig(PCA_FIGURE_OUT, dpi=300)
    plt.close(fig)


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No data._"
    text_df = df.copy()
    for column in text_df.columns:
        text_df[column] = text_df[column].map(lambda value: str(value).replace("|", "\\|"))
    header = "| " + " | ".join(text_df.columns) + " |"
    separator = "| " + " | ".join("---" for _ in text_df.columns) + " |"
    rows = [
        "| " + " | ".join(row[column] for column in text_df.columns) + " |"
        for _, row in text_df.iterrows()
    ]
    return "\n".join([header, separator] + rows)


def candidate_diagnostic_table(sample_qc: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "sample_id",
        "source_cohort",
        "response_pcr_vs_rd",
        "expression_mean",
        "expression_median",
        "expression_sd",
        "expression_iqr",
        "robust_z_expression_iqr",
        "PC1",
        "robust_z_PC1",
        "PC2",
        "robust_z_PC2",
        "median_correlation_all_samples",
        "median_correlation_within_cohort",
        "robust_z_median_correlation_within_cohort",
        "distribution_flag",
        "pca_extreme_flag",
    ]
    return sample_qc.loc[sample_qc["qc_candidate_flag"] == "True", columns].copy()


def write_report(
    status: str,
    checksums: dict[str, str],
    metadata: pd.DataFrame,
    cohort_table: pd.DataFrame,
    age_summary: dict[str, object],
    missing_rows: list[dict[str, object]],
    sample_qc: pd.DataFrame,
    explained: list[float],
    incidents: list[str],
) -> None:
    response_counts = Counter(metadata["response_pcr_vs_rd"])
    primary = metadata.loc[metadata["primary_modeling_eligible"] == "True"]
    primary_counts = Counter(primary["response_pcr_vs_rd"])
    qc_candidates = sample_qc.loc[sample_qc["qc_candidate_flag"] == "True", "sample_id"].tolist()
    previous_set = set(PREVIOUS_QC_CANDIDATES)
    current_set = set(qc_candidates)
    distribution_ids = sample_qc.loc[sample_qc["distribution_flag"] == "True", "sample_id"].tolist()
    pca_ids = sample_qc.loc[sample_qc["pca_extreme_flag"] == "True", "sample_id"].tolist()
    duplicate_ids = sample_qc.loc[sample_qc["exact_duplicate_profile_flag"] == "True", "sample_id"].tolist()
    low_corr_ids = sample_qc.loc[
        pd.to_numeric(sample_qc["robust_z_median_correlation_within_cohort"], errors="coerce").lt(-5), "sample_id"
    ].tolist()
    expr_stats = sample_qc[["expression_mean", "expression_median", "expression_sd", "expression_iqr"]].apply(pd.to_numeric)

    lines = [
        "# GSE25066: EDA y control de calidad",
        "",
        "## 1. Estado",
        "",
        f"**{status}**",
        "",
        "## 2. Inputs y checksums",
        "",
    ]
    lines += [f"- `{key}`: `{value}`" for key, value in checksums.items()]
    lines += [
        "",
        "## 3. Poblaciones utilizadas",
        "",
        f"- Endpoint-known: {len(metadata)} muestras.",
        f"- Primaria para PCA/QC de expresion: {len(primary)} muestras.",
        "",
        "## 4. Balance de respuesta y prevalencia",
        "",
        f"- Endpoint-known: pCR={response_counts['pCR']}, RD={response_counts['RD']}, prevalencia pCR={100 * response_counts['pCR'] / len(metadata):.2f}%.",
        f"- Primaria: pCR={primary_counts['pCR']}, RD={primary_counts['RD']}, prevalencia pCR={100 * primary_counts['pCR'] / len(primary):.2f}%.",
        "",
        "## 5. Distribucion por cohortes",
        "",
        markdown_table(cohort_table),
        "",
        "## 6. Resumen de ER, PR, HER2 y edad",
        "",
        "- ER endpoint-known: " + ", ".join(f"{k}={v}" for k, v in sorted(normalized_counter(metadata["er_status"]).items())),
        "- PR endpoint-known: " + ", ".join(f"{k}={v}" for k, v in sorted(normalized_counter(metadata["pr_status"]).items())),
        "- HER2 endpoint-known: " + ", ".join(f"{k}={v}" for k, v in sorted(normalized_counter(metadata["her2_status"]).items())),
        f"- Edad: min={age_summary['min']:.2f}, Q1={age_summary['q1']:.2f}, mediana={age_summary['median']:.2f}, Q3={age_summary['q3']:.2f}, max={age_summary['max']:.2f}.",
        f"- Edad missing={age_summary['missing']}, no numerica={age_summary['non_numeric']}, fuera de 18-100={len(age_summary['outside_18_100_ids'])}.",
        "- RCB, DRFS evento y DRFS tiempo estan disponibles como inventario, pero no se han usado como predictores ni para excluir muestras.",
        "",
        "## 7. Missingness de variables clinicas",
        "",
        markdown_table(pd.DataFrame(missing_rows)),
        "",
        "## 8. Resumen global de expresion",
        "",
    ]
    for column in expr_stats.columns:
        lines.append(
            f"- {column}: min={expr_stats[column].min():.4f}, Q1={expr_stats[column].quantile(0.25):.4f}, "
            f"mediana={expr_stats[column].median():.4f}, Q3={expr_stats[column].quantile(0.75):.4f}, max={expr_stats[column].max():.4f}."
        )
    lines += [
        "",
        "No se interpreta ninguna sonda concreta en este análisis.",
        "",
        "## 9. PCA",
        "",
        "- PCA descriptiva sobre 482 muestras primarias y 22.283 sondas.",
        "- Implementacion estandar: `sklearn.decomposition.PCA(n_components=10, svd_solver=\"randomized\", random_state=20260608)`.",
        "- PCA realiza el centrado interno; no se escalo a varianza unitaria.",
        "- La PCA no es transformacion del pipeline predictivo y no podra reutilizarse como seleccion o reduccion dimensional en el modelado.",
        "- Varianza explicada PC1-PC10: "
        + ", ".join(f"PC{i + 1}={value:.2f}%" for i, value in enumerate(explained)),
        "- La estructura visible se documenta de forma neutral por cohorte, ER y respuesta en la figura correspondiente, sin afirmar separabilidad predictiva.",
        "",
        "## 10. Control de muestras",
        "",
        f"- Candidatos QC previos: {', '.join(PREVIOUS_QC_CANDIDATES)}.",
        f"- Candidatos QC: {', '.join(qc_candidates) if qc_candidates else 'ninguno'}.",
        f"- Nuevos candidatos respecto a la ejecucion previa: {', '.join(sorted(current_set - previous_set)) if current_set - previous_set else 'ninguno'}.",
        f"- Candidatos previos no retenidos: {', '.join(sorted(previous_set - current_set)) if previous_set - current_set else 'ninguno'}.",
        f"- distribution_flag: {len(distribution_ids)}.",
        f"- pca_extreme_flag: {len(pca_ids)}.",
        f"- perfiles duplicados exactos: {len(duplicate_ids)}.",
        f"- correlacion mediana dentro de cohorte con robust z < -5: {len(low_corr_ids)}.",
        f"- candidatos QC totales: {len(qc_candidates)}.",
        "",
        "### Revision tecnica focalizada de candidatos QC",
        "",
        markdown_table(candidate_diagnostic_table(sample_qc)),
        "",
        "- Los candidatos se conservan provisionalmente.",
        "- Ninguna muestra se excluye por PCA o dispersion global.",
        "- Cualquier exclusion requeriria evidencia tecnica adicional preespecificada.",
        "- Podra realizarse una sensibilidad posterior con y sin candidatos, sin escoger post hoc el resultado mas favorable.",
        "",
        "## 11. Variables deliberadamente no utilizadas",
        "",
        "- RCB.",
        "- DRFS.",
        "- Anotacion genica para seleccion.",
        "- Respuesta para seleccionar sondas.",
        "- Correlaciones tecnicas entre perfiles: solo diagnostico QC, no predictor.",
        "",
        "## 12. Outputs generados",
        "",
        f"- `{rel(REPORT_OUT)}`",
        f"- `{rel(COHORT_SUMMARY_OUT)}`",
        f"- `{rel(SAMPLE_QC_OUT)}`",
        f"- `{rel(CLINICAL_FIGURE_OUT)}`",
        f"- `{rel(PCA_FIGURE_OUT)}`",
        "",
        "## 13. Analisis no realizados",
        "",
        "- Filtrado de sondas para modelado.",
        "- Escalado para modelado.",
        "- Seleccion de features.",
        "- Entrenamiento.",
        "- Validacion cruzada.",
        "- Analisis diferencial.",
        "- Enriquecimiento.",
        "",
        "## Incidencias",
        "",
    ]
    lines += [f"- {incident}" for incident in incidents] if incidents else ["- Ninguna incidencia estructural."]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dependency_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit-learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "requests": requests.__version__,
    }



def main() -> int:
    print("Iniciando EDA y control de calidad para GSE25066.")
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    checksums = {
        rel(METADATA_PATH): sha256_file(METADATA_PATH),
        rel(EXPRESSION_PATH): sha256_file(EXPRESSION_PATH),
    }
    metadata, expression, probe_ids, profile_hashes = load_inputs()
    age, age_summary = numeric_age(metadata)
    missing_rows = missingness_table(metadata)
    cohort_table = cohort_summary(metadata)
    sample_qc, explained, incidents = run_qc_and_pca(metadata, expression, probe_ids, profile_hashes)

    status_issues: list[str] = []
    missing_lookup = {row["variable"]: row for row in missing_rows}
    for variable, expected in EXPECTED_MISSINGNESS.items():
        row = missing_lookup[variable]
        if row["n_missing"] != expected["missing"] or row["n_indeterminate"] != expected["indeterminate"]:
            status_issues.append(
                f"{variable} missing/indeterminate mismatch: observed {row['n_missing']}/{row['n_indeterminate']}, expected {expected['missing']}/{expected['indeterminate']}"
            )

    duplicate_n = int((sample_qc["exact_duplicate_profile_flag"] == "True").sum())
    if duplicate_n:
        status_issues.append(f"{duplicate_n} samples with exact duplicate profiles")
    if age_summary["non_numeric"]:
        status_issues.append(f"{age_summary['non_numeric']} non-numeric age values")
    if age_summary["non_finite"]:
        status_issues.append(f"{age_summary['non_finite']} non-finite age values")
    if age_summary["outside_18_100_ids"]:
        status_issues.append("age values outside 18-100: " + ", ".join(age_summary["outside_18_100_ids"]))
    low_corr_ids = sample_qc.loc[
        pd.to_numeric(sample_qc["robust_z_median_correlation_within_cohort"], errors="coerce").lt(-5), "sample_id"
    ].tolist()
    if low_corr_ids:
        status_issues.append("samples with robust_z_median_correlation_within_cohort < -5: " + ", ".join(low_corr_ids))

    status = "REVIEW_REQUIRED" if status_issues else "ANALYSIS_COMPLETE"
    incidents.extend(status_issues)
    if not status_issues:
        incidents.append("Candidatos QC documentados y retenidos; no hay incidencias estructurales.")

    cohort_table.to_csv(COHORT_SUMMARY_OUT, sep="\t", index=False)
    sample_qc.to_csv(SAMPLE_QC_OUT, sep="\t", index=False)
    write_clinical_figure(metadata, age, missing_rows)
    write_pca_figure(sample_qc, explained)
    write_report(status, checksums, metadata, cohort_table, age_summary, missing_rows, sample_qc, explained, incidents)

    distribution_n = int((sample_qc["distribution_flag"] == "True").sum())
    pca_n = int((sample_qc["pca_extreme_flag"] == "True").sum())
    candidate_ids = sample_qc.loc[sample_qc["qc_candidate_flag"] == "True", "sample_id"].tolist()
    print("Versions: " + ", ".join(f"{key}={value}" for key, value in dependency_versions().items()))
    print(f"Dimensions: metadata=488 rows; expression=488 x 22283 probes; primary=482")
    for variable in ["er_status", "pr_status", "her2_status", "grade"]:
        row = missing_lookup[variable]
        print(f"{variable}: missing={row['n_missing']} ({row['pct_missing']:.2f}%), indeterminate={row['n_indeterminate']}")
    print(f"PCA variance: PC1={explained[0]:.2f}%, PC2={explained[1]:.2f}%")
    print("Previous QC candidates: " + ", ".join(PREVIOUS_QC_CANDIDATES))
    print("Current QC candidates: " + (", ".join(candidate_ids) if candidate_ids else "none"))
    diagnostic = candidate_diagnostic_table(sample_qc)
    if diagnostic.empty:
        print("Candidate correlations: none")
    else:
        print("Candidate correlations:")
        for _, row in diagnostic.iterrows():
            print(
                f"- {row['sample_id']}: median_all={row['median_correlation_all_samples']}, "
                f"median_within_cohort={row['median_correlation_within_cohort']}, "
                f"robust_z_within={row['robust_z_median_correlation_within_cohort']}"
            )
    print(f"Distribution flags: {distribution_n}")
    print(f"PCA flags: {pca_n}")
    print(f"Exact duplicate profiles: {duplicate_n}")
    print("Exclusions applied: 0")
    print(f"Final state: {status}")
    print(f"STATUS: {status}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
