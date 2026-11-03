"""Analyte-level (column) quality control and the detectability filter."""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

MAD_CONST = 1.4826


def limits_of_detection(buffer_rfu: pd.DataFrame, mad_multiplier: float = 4.9) -> pd.Series:
    """Per-analyte LoD estimated from buffer (blank) wells.

    LoD = median(buffer) + k * MAD(buffer), the blank-based estimator used for
    SomaScan detectability work (k = 4.9 corresponds to the SomaLogic convention).
    """
    med = buffer_rfu.median(axis=0)
    mad = (buffer_rfu - med).abs().median(axis=0) * MAD_CONST
    return med + mad_multiplier * mad


def detection_matrix(rfu: pd.DataFrame, lod: pd.Series) -> pd.DataFrame:
    return rfu.gt(lod.reindex(rfu.columns), axis=1)


def replicate_cv(rfu: pd.DataFrame, plate: pd.Series | None) -> pd.Series:
    """Median intra-plate coefficient of variation from pooled QC replicates."""
    if rfu.shape[0] < 3:
        return pd.Series(np.nan, index=rfu.columns)
    if plate is None or plate.nunique() <= 1:
        cv = rfu.std(axis=0, ddof=1) / rfu.mean(axis=0)
        return cv
    per_plate = []
    for _, idx in plate.groupby(plate).groups.items():
        sub = rfu.loc[list(idx)]
        if sub.shape[0] >= 3:
            per_plate.append(sub.std(axis=0, ddof=1) / sub.mean(axis=0))
    if not per_plate:
        return pd.Series(np.nan, index=rfu.columns)
    return pd.concat(per_plate, axis=1).median(axis=1)


def analyte_qc(analytes: pd.DataFrame, study_rfu: pd.DataFrame, groups: pd.Series,
               lod: pd.Series | None, qc_rfu: pd.DataFrame | None,
               qc_plate: pd.Series | None, params: Dict) -> Tuple[pd.DataFrame, pd.DataFrame | None]:
    """Annotate every analyte with QC flags and the filtering decision.

    Returns the per-analyte table and the per-group detection-rate table.
    """
    p = params
    a = pd.DataFrame(index=analytes.index)
    keep_cols = [c for c in ("SomaId", "Target", "TargetFullName", "UniProt",
                             "EntrezGeneSymbol", "EntrezGeneID", "Organism", "Type",
                             "Dilution", "Units", "ColCheck") if c in analytes.columns]
    for c in keep_cols:
        a[c] = analytes[c]

    # --- content class -------------------------------------------------
    atype = analytes.get("Type", pd.Series("Protein", index=analytes.index)).astype(str)
    a["flag_not_protein"] = ~atype.isin(p["keep_analyte_types"])
    org = analytes.get("Organism", pd.Series("Human", index=analytes.index)).astype(str)
    a["flag_organism"] = (~org.isin(p["keep_organisms"])) if p["keep_organisms"] else False

    # --- SomaLogic column QC -------------------------------------------
    colcheck = analytes.get("ColCheck", pd.Series("PASS", index=analytes.index))
    a["flag_colcheck"] = colcheck.astype(str).str.upper().eq("FLAG")

    ratio_cols = [c for c in analytes.columns if c.startswith("CalQcRatio")]
    if ratio_cols:
        ratios = analytes[ratio_cols].apply(pd.to_numeric, errors="coerce")
        out = (ratios < p["cal_qc_ratio_min"]) | (ratios > p["cal_qc_ratio_max"])
        a["n_plates_calqc_out"] = out.sum(axis=1)
        a["flag_cal_qc_ratio"] = out.any(axis=1)
    else:
        a["n_plates_calqc_out"] = 0
        a["flag_cal_qc_ratio"] = False

    # --- QC-replicate CV ------------------------------------------------
    if qc_rfu is not None and qc_rfu.shape[0] >= 3:
        cv = replicate_cv(qc_rfu, qc_plate)
        a["qc_replicate_cv"] = cv.reindex(a.index)
        a["flag_high_cv"] = a["qc_replicate_cv"] > p["max_replicate_cv"]
    else:
        a["qc_replicate_cv"] = np.nan
        a["flag_high_cv"] = False

    # --- detectability ---------------------------------------------------
    det_table = None
    if lod is not None and p["detectability_filter"]:
        a["lod_rfu"] = lod.reindex(a.index)
        det = detection_matrix(study_rfu, lod)
        rates = det.groupby(groups.reindex(det.index).values, observed=True).mean()
        det_table = rates.T
        det_table.columns = [f"detect_rate::{c}" for c in det_table.columns]
        a["detect_rate_overall"] = det.mean(axis=0).reindex(a.index)
        passes = (det_table >= p["min_detect_frac"])
        a["n_groups_detected"] = passes.sum(axis=1).reindex(a.index)
        if p["detect_rule"] == "all":
            keep_det = passes.all(axis=1)
        else:
            keep_det = passes.any(axis=1)
        a["flag_low_detection"] = ~keep_det.reindex(a.index).fillna(False)
        for c in det_table.columns:
            a[c] = det_table[c]
    else:
        a["lod_rfu"] = np.nan
        a["detect_rate_overall"] = np.nan
        a["n_groups_detected"] = np.nan
        a["flag_low_detection"] = False

    drop_on: List[str] = [c for c in p["analyte_exclude_on"] if c in a.columns]
    a["flag_any"] = a[[c for c in a.columns if c.startswith("flag_")]].any(axis=1)
    a["exclude"] = a[drop_on].any(axis=1) if drop_on else pd.Series(False, index=a.index)
    a["exclude_reason"] = [
        ";".join(c.replace("flag_", "") for c in drop_on if bool(a.at[i, c])) for i in a.index
    ]
    return a, det_table
