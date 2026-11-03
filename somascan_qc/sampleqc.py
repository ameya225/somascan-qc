"""Sample-level (row) quality control.

Criteria follow SomaLogic's published acceptance rules for ANML-normalized data,
plus two standard data-driven outlier checks:

  * ``RowCheck``            - SomaLogic per-sample normalization acceptance flag.
  * normalization scale     - median scale factor per dilution group must sit
                              inside the acceptance window (default 0.4-2.5);
                              same window applied to the hybridization control
                              normalization scale factor.
  * ``ANMLFractionUsed``    - fraction of analytes inside the reference
                              distribution used by ANML; a low value means the
                              sample deviates strongly from the reference.
  * detection rate          - fraction of analytes above the buffer-derived LoD.
  * median-correlation      - Spearman-free: Pearson correlation of the sample's
                              log2 profile to the study median profile; low
                              correlation marks global outliers.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    return df[cols].apply(pd.to_numeric, errors="coerce") if cols else pd.DataFrame(index=df.index)


def sample_qc(samples: pd.DataFrame, rfu: pd.DataFrame, params: dict,
              detect_rate: pd.Series | None = None) -> pd.DataFrame:
    """Return a per-sample QC table with one boolean column per criterion."""
    p = params
    qc = pd.DataFrame(index=samples.index)
    qc["SampleType"] = samples.get("SampleType", pd.Series("Sample", index=samples.index))
    if "PlateId" in samples:
        qc["PlateId"] = samples["PlateId"]

    ns_cols = [c for c in samples.columns if c.startswith("NormScale")]
    hyb_cols = [c for c in samples.columns if c.startswith("HybControlNormScale")]
    anml_cols = [c for c in samples.columns if c.startswith("ANMLFractionUsed")]

    ns = _numeric(samples, ns_cols)
    hyb = _numeric(samples, hyb_cols)
    anml = _numeric(samples, anml_cols)
    for frame in (ns, hyb, anml):
        for c in frame.columns:
            qc[c] = frame[c]

    lo, hi = p["norm_scale_min"], p["norm_scale_max"]
    qc["flag_rowcheck"] = (samples.get("RowCheck", pd.Series("PASS", index=samples.index))
                           .astype(str).str.upper().eq("FLAG"))
    qc["flag_norm_scale"] = (((ns < lo) | (ns > hi)).any(axis=1) if not ns.empty
                             else pd.Series(False, index=samples.index))
    qc["flag_hyb_scale"] = (((hyb < p["hyb_scale_min"]) | (hyb > p["hyb_scale_max"])).any(axis=1)
                            if not hyb.empty else pd.Series(False, index=samples.index))
    qc["flag_anml_fraction"] = ((anml < p["anml_fraction_min"]).any(axis=1) if not anml.empty
                                else pd.Series(False, index=samples.index))

    log2 = np.log2(rfu.clip(lower=1e-6))
    qc["median_log2_rfu"] = log2.median(axis=1)
    is_study = qc["SampleType"].eq("Sample")
    ref = log2.loc[is_study].median(axis=0) if is_study.any() else log2.median(axis=0)
    X = log2.values - log2.values.mean(axis=1, keepdims=True)
    refc = ref.values - ref.values.mean()
    denom = np.sqrt((X ** 2).sum(axis=1) * (refc ** 2).sum())
    qc["cor_to_median"] = np.where(denom > 0, (X @ refc) / np.where(denom > 0, denom, 1), np.nan)
    qc["flag_low_correlation"] = qc["cor_to_median"] < p["min_correlation_to_median"]

    if detect_rate is not None:
        qc["detect_rate"] = detect_rate.reindex(qc.index)
        qc["flag_low_detection"] = qc["detect_rate"] < p["sample_min_detect_rate"]
    else:
        qc["detect_rate"] = np.nan
        qc["flag_low_detection"] = False

    # Robust intensity outlier: MAD-based z of the per-sample median, study samples only.
    med = qc.loc[is_study, "median_log2_rfu"]
    centre = med.median()
    mad = (med - centre).abs().median() * 1.4826
    z = (qc["median_log2_rfu"] - centre) / mad if mad > 0 else pd.Series(0.0, index=qc.index)
    qc["median_intensity_z"] = z
    qc["flag_intensity_outlier"] = z.abs() > p["max_intensity_mad_z"]

    notes = samples.get("SampleNotes", pd.Series("", index=samples.index)).astype(str)
    qc["SampleNotes"] = notes
    qc["flag_sample_note"] = notes.str.strip().ne("")

    active = [c for c in p["sample_exclude_on"] if c in qc.columns]
    qc["flag_any"] = qc[[c for c in qc.columns if c.startswith("flag_")]].any(axis=1)
    qc["exclude"] = qc[active].any(axis=1) if active else pd.Series(False, index=qc.index)
    qc.loc[~is_study, "exclude"] = False  # control wells are handled by sample-type selection
    reasons = []
    for i in qc.index:
        r = [c.replace("flag_", "") for c in active if bool(qc.at[i, c])]
        reasons.append(";".join(r))
    qc["exclude_reason"] = reasons
    return qc
