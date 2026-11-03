"""Parameter defaults, validation and provenance hashing."""
from __future__ import annotations

import hashlib
import json
from typing import Dict

DEFAULTS: Dict = {
    # ---- input / output -------------------------------------------------
    "input_adat": None,
    "output_dir": "somascan_qc_out",
    "sample_types_to_keep": ["Sample"],
    "group_col": None,              # any sample column in the ADAT; None = one pooled group
    "output_prefix": "somascan",

    # ---- sample QC ------------------------------------------------------
    "norm_scale_min": 0.4,
    "norm_scale_max": 2.5,
    "hyb_scale_min": 0.4,
    "hyb_scale_max": 2.5,
    "anml_fraction_min": 0.3,
    "min_correlation_to_median": 0.8,
    "max_intensity_mad_z": 5.0,
    "sample_min_detect_rate": 0.5,
    "sample_exclude_on": ["flag_rowcheck", "flag_norm_scale", "flag_hyb_scale",
                          "flag_anml_fraction", "flag_low_correlation",
                          "flag_intensity_outlier", "flag_low_detection"],

    # ---- analyte QC -----------------------------------------------------
    "keep_analyte_types": ["Protein"],
    "keep_organisms": ["Human"],      # [] disables the organism filter
    "cal_qc_ratio_min": 0.8,
    "cal_qc_ratio_max": 1.2,
    "max_replicate_cv": 0.2,
    "detectability_filter": True,
    "lod_mad_multiplier": 4.9,
    "min_detect_frac": 0.5,
    "detect_rule": "any",             # "any" | "all" groups must pass
    "analyte_exclude_on": ["flag_not_protein", "flag_organism", "flag_low_detection"],

    # ---- transform / output --------------------------------------------
    "log_transform": "log2",          # "log2" | "log10" | "none"
    "write_rfu_matrix": True,
    "compress_outputs": False,       # write .tsv.gz instead of .tsv
    "matrix_column_id": "sample_key",  # sample column used as matrix column names
    "plate_median_center": False,
    "report_filename": "qc_report.html",
    "deterministic_report": False,   # omit timestamp/runtime so the HTML is byte-reproducible
}

_NUMERIC = {k: (0.0, 1e6) for k in ("norm_scale_min", "norm_scale_max", "hyb_scale_min",
                                    "hyb_scale_max", "lod_mad_multiplier", "max_intensity_mad_z")}
_FRACTION = ["anml_fraction_min", "min_correlation_to_median", "sample_min_detect_rate",
             "min_detect_frac", "max_replicate_cv", "cal_qc_ratio_min", "cal_qc_ratio_max"]


def merge_params(user: Dict | None) -> Dict:
    """Merge user settings over the defaults and validate."""
    p = dict(DEFAULTS)
    unknown = sorted(set(user or {}) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"unknown parameter(s) in config: {', '.join(unknown)}. "
                         f"Valid keys: {', '.join(sorted(DEFAULTS))}")
    p.update({k: v for k, v in (user or {}).items() if v is not None})
    for k in _FRACTION:
        if not (0 <= float(p[k]) <= 10):
            raise ValueError(f"{k} must be >= 0 (got {p[k]})")
    if p["norm_scale_min"] >= p["norm_scale_max"]:
        raise ValueError("norm_scale_min must be < norm_scale_max")
    if p["detect_rule"] not in ("any", "all"):
        raise ValueError("detect_rule must be 'any' or 'all'")
    if p["log_transform"] not in ("log2", "log10", "none"):
        raise ValueError("log_transform must be log2, log10 or none")
    for k in ("sample_exclude_on", "analyte_exclude_on", "keep_analyte_types",
              "keep_organisms", "sample_types_to_keep"):
        if not isinstance(p[k], list):
            raise ValueError(f"{k} must be a list")
    return p


def param_hash(p: Dict) -> str:
    """Stable digest of the effective parameters (input path excluded)."""
    q = {k: v for k, v in p.items() if k not in ("input_adat", "output_dir")}
    return hashlib.sha256(json.dumps(q, sort_keys=True, default=str).encode()).hexdigest()[:12]
