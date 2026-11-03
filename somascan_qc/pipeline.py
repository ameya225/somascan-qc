"""End-to-end QC pipeline: ADAT in, filtered protein x sample matrix + report out."""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import time
from typing import Dict, List

import numpy as np
import pandas as pd

from . import __version__, report
from .adat import Adat, read_adat
from .analyteqc import analyte_qc, limits_of_detection
from .params import merge_params, param_hash
from .sampleqc import sample_qc


def _sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _write(df: pd.DataFrame, path: str, index_label: str, compress: bool = False) -> str:
    if compress:
        path = path + ".gz"
    df.to_csv(path, sep="\t", index_label=index_label, float_format="%.6g",
              compression="gzip" if compress else None)
    return path


def resolve_groups(adat: Adat, samples: pd.DataFrame, group_col: str | None,
                   warnings: List[str]) -> pd.Series:
    """Resolve the grouping variable from a sample-annotation column of the ADAT."""
    if group_col in (None, "", "none"):
        warnings.append(
            "No <code>group_col</code> set: detectability was evaluated on all retained "
            "samples as a single pooled group. Set <code>group_col</code> to any populated "
            "ADAT sample column to filter per group.")
        return pd.Series("all_samples", index=samples.index, name="group")
    if group_col not in samples.columns:
        cands = ", ".join(adat.grouping_candidates()) or "none found"
        raise ValueError(f"group_col '{group_col}' is not a sample column in the ADAT. "
                         f"Usable columns: {cands}")
    g = samples[group_col].astype(str).str.strip().replace("", "unassigned")
    if g.nunique() < 2:
        warnings.append(f"<code>group_col='{group_col}'</code> has a single level - "
                        "detectability filtering is effectively ungrouped.")
    return g.rename("group")


def run(params: Dict) -> Dict:
    """Execute the pipeline. Returns the run manifest."""
    t0 = time.time()
    p = merge_params(params)
    if not p["input_adat"]:
        raise ValueError("input_adat is required")
    outdir = p["output_dir"]
    os.makedirs(outdir, exist_ok=True)
    warnings: List[str] = []
    timings: Dict[str, float] = {}

    # 1 -------------------------------------------------------------- read
    t = time.time()
    adat = read_adat(p["input_adat"])
    timings["read_adat_s"] = time.time() - t

    stype = adat.samples.get("SampleType", pd.Series("Sample", index=adat.samples.index))
    is_study = stype.isin(p["sample_types_to_keep"])
    if not is_study.any():
        raise ValueError(f"no samples of type {p['sample_types_to_keep']} in {p['input_adat']}")
    is_buffer = stype.eq("Buffer")
    is_qc = stype.eq("QC")

    # 2 ------------------------------------------------- limits of detection
    t = time.time()
    lod = None
    if p["detectability_filter"]:
        if is_buffer.sum() >= 3:
            lod = limits_of_detection(adat.rfu.loc[is_buffer.values],
                                      p["lod_mad_multiplier"])
        else:
            warnings.append(
                f"Only {int(is_buffer.sum())} buffer well(s) in this ADAT; a blank-based limit "
                "of detection could not be estimated, so the detectability filter was skipped. "
                "Analytes were still flagged on all other criteria.")
            p = dict(p, detectability_filter=False)
    timings["lod_s"] = time.time() - t

    # 3 ------------------------------------------------------- sample QC
    t = time.time()
    det_rate = None
    if lod is not None:
        det_rate = adat.rfu.gt(lod.reindex(adat.rfu.columns), axis=1).mean(axis=1)
    sqc = sample_qc(adat.samples, adat.rfu, p, det_rate)
    keep_samples = is_study.values & (~sqc["exclude"].values)
    timings["sample_qc_s"] = time.time() - t

    if keep_samples.sum() == 0:
        raise ValueError(
            "every study sample was excluded by sample QC - inspect sample_qc_all.tsv and "
            "loosen sample_exclude_on thresholds")
    study = adat.samples.loc[keep_samples]
    study_rfu = adat.rfu.loc[keep_samples]
    groups = resolve_groups(adat, study, p["group_col"], warnings)

    # 4 ------------------------------------------------------ analyte QC
    t = time.time()
    qc_rfu = adat.rfu.loc[is_qc.values] if is_qc.sum() >= 3 else None
    qc_plate = adat.samples.loc[is_qc.values, "PlateId"] if (
        qc_rfu is not None and "PlateId" in adat.samples) else None
    aqc, det_table = analyte_qc(adat.analytes, study_rfu, groups, lod, qc_rfu, qc_plate, p)
    keep_analytes = ~aqc["exclude"].values
    if keep_analytes.sum() == 0:
        raise ValueError("every analyte was excluded - loosen analyte_exclude_on / "
                         "min_detect_frac, or check that the ADAT is ANML-normalized")
    timings["analyte_qc_s"] = time.time() - t

    # 5 ------------------------------------------- matrix build + transform
    t = time.time()
    mat_rfu = study_rfu.loc[:, keep_analytes].T          # proteins x samples
    mat_rfu.index.name = "SeqId"
    col_id = p["matrix_column_id"]
    if col_id not in ("sample_key", None) :
        if col_id not in study.columns:
            raise ValueError(f"matrix_column_id '{col_id}' is not a sample column in the ADAT")
        new = study[col_id].astype(str)
        if not new.is_unique:
            raise ValueError(f"matrix_column_id '{col_id}' is not unique across retained samples")
        mat_rfu.columns = new.values
    if p["log_transform"] == "log2":
        mat = np.log2(mat_rfu.clip(lower=1e-6))
    elif p["log_transform"] == "log10":
        mat = np.log10(mat_rfu.clip(lower=1e-6))
    else:
        mat = mat_rfu.copy()
    if p["plate_median_center"] and "PlateId" in study:
        plate = study["PlateId"].astype(str)
        grand = mat.median(axis=1)
        for pl in plate.unique():
            cols = list(plate.index[plate.eq(pl)])
            mat[cols] = mat[cols].sub(mat[cols].median(axis=1), axis=0).add(grand, axis=0)
        warnings.append("Per-plate median centering was applied to the output matrix "
                        "(each analyte's plate median was shifted onto its overall median).")
    timings["matrix_s"] = time.time() - t

    # 6 ------------------------------------------------------------ write
    t = time.time()
    pre = os.path.join(outdir, p["output_prefix"])
    gz = bool(p["compress_outputs"])
    suffix = {"log2": "log2", "log10": "log10", "none": "rfu"}[p["log_transform"]]
    outputs: Dict[str, str] = {}
    outputs["matrix"] = _write(mat, f"{pre}_protein_x_sample_{suffix}.tsv", "SeqId", gz)
    if p["write_rfu_matrix"] and p["log_transform"] != "none":
        outputs["matrix_rfu"] = _write(mat_rfu, f"{pre}_protein_x_sample_rfu.tsv", "SeqId", gz)
    feat = aqc.loc[keep_analytes]
    outputs["feature_metadata"] = _write(feat, f"{pre}_feature_metadata.tsv", "SeqId", gz)
    outputs["analyte_qc_all"] = _write(aqc, f"{pre}_analyte_qc_all.tsv", "SeqId", gz)
    smeta_out = study.copy()
    smeta_out["group"] = groups
    smeta_out["matrix_column"] = mat.columns
    for c in ("median_log2_rfu", "cor_to_median", "detect_rate", "median_intensity_z"):
        smeta_out[c] = sqc.loc[study.index, c]
    outputs["sample_metadata"] = _write(smeta_out, f"{pre}_sample_metadata.tsv", "sample_key", gz)
    outputs["sample_qc_all"] = _write(sqc, f"{pre}_sample_qc_all.tsv", "sample_key", gz)
    timings["write_s"] = time.time() - t

    # 7 ----------------------------------------------------------- report
    t = time.time()
    ctx = _report_context(adat, p, sqc, aqc, det_table, groups, mat, mat_rfu, lod,
                          is_study, is_buffer, is_qc, keep_samples, keep_analytes,
                          warnings, time.time() - t0)
    html = report.build_report(ctx)
    f = os.path.join(outdir, p["report_filename"])
    with open(f, "w", encoding="utf-8") as fh:
        fh.write(html)
    outputs["report"] = f
    timings["report_s"] = time.time() - t

    # 8 --------------------------------------------------------- manifest
    elapsed = time.time() - t0
    manifest = {
        "pipeline": "somascan_qc", "version": __version__,
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        "input": {"path": os.path.abspath(p["input_adat"]), "sha256": adat.sha256,
                  "n_samples": adat.n_samples, "n_analytes": adat.n_analytes,
                  "sample_types": adat.sample_types(),
                  "adat_header": {k: adat.header.get(k, "") for k in
                                  ("AssayVersion", "ProcessSteps", "NormalizationAlgorithm",
                                   "StudyMatrix", "CreatedDate", "AdatId")}},
        "parameters": p, "param_hash": param_hash(p),
        "results": {
            "samples_in": int(adat.n_samples),
            "study_samples": int(is_study.sum()),
            "samples_excluded": int(is_study.sum() - keep_samples.sum()),
            "samples_out": int(keep_samples.sum()),
            "analytes_in": int(adat.n_analytes),
            "analytes_out": int(keep_analytes.sum()),
            "groups": {str(k): int(v) for k, v in groups.value_counts().items()},
        },
        "timings_s": {k: round(v, 3) for k, v in timings.items()} | {"total_s": round(elapsed, 3)},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "pandas": pd.__version__},
        "outputs": {k: {"path": os.path.abspath(v), "sha256": _sha(v)} for k, v in outputs.items()},
        "warnings": [w.replace("<code>", "").replace("</code>", "") for w in warnings],
    }
    mpath = os.path.join(outdir, f"{p['output_prefix']}_run_manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    manifest["manifest_path"] = mpath
    return manifest


def _pct(n, d):
    return f"{100.0 * n / d:.1f}" if d else "0.0"


def _report_context(adat, p, sqc, aqc, det_table, groups, mat, mat_rfu, lod,
                    is_study, is_buffer, is_qc, keep_samples, keep_analytes,
                    warnings, elapsed) -> Dict:
    report._style()
    study_sqc = sqc[sqc.SampleType.eq("Sample")]
    n_study = int(is_study.sum())
    n_keep = int(keep_samples.sum())

    figs = {
        "norm_scale": report.fig_norm_scale(sqc, p),
        "anml": report.fig_anml(sqc, p),
        "sample_intensity": report.fig_sample_intensity(sqc),
        "sample_hists": report.fig_sample_hists(sqc, p),
        "intensity_by_type": report.fig_intensity_by_type(
            adat.rfu, adat.samples.get("SampleType", pd.Series("Sample", index=adat.samples.index))),
        "analyte_panels": report.fig_analyte_panels(aqc, p),
    }
    log2_mat = mat.T if p["log_transform"] == "log2" else np.log2(mat_rfu.T.clip(lower=1e-6))
    plate = adat.samples.loc[keep_samples, "PlateId"] if "PlateId" in adat.samples else \
        pd.Series("all", index=log2_mat.index)
    figs["pca"] = report.fig_pca(log2_mat, plate, groups)
    figs["plate_effect"] = report.fig_plate_effect(log2_mat, plate)

    steps = [("all analytes in ADAT", int(adat.n_analytes))]
    cum = ~aqc["flag_not_protein"].values
    steps.append((f"type in {p['keep_analyte_types']}", int(cum.sum()))) if True else None
    if p["keep_organisms"]:
        cum = cum & ~aqc["flag_organism"].values
        steps.append((f"organism in {p['keep_organisms']}", int(cum.sum())))
    if p["detectability_filter"]:
        cum2 = cum & ~aqc["flag_low_detection"].values
        steps.append((f"detected in >={p['min_detect_frac']:.0%} of samples "
                      f"({p['detect_rule']} group)", int(cum2.sum())))
    steps.append(("retained in output matrix", int(keep_analytes.sum())))
    figs["waterfall"] = report.fig_waterfall(steps)

    sflags = [c for c in sqc.columns if c.startswith("flag_") and c != "flag_any"]
    sft = pd.DataFrame({
        "criterion": [c.replace("flag_", "") for c in sflags],
        "samples flagged": [int(study_sqc[c].sum()) for c in sflags],
        "% of study samples": [_pct(int(study_sqc[c].sum()), n_study) for c in sflags],
        "used for exclusion": ["yes" if c in p["sample_exclude_on"] else "flag only"
                               for c in sflags]})
    aflags = [c for c in aqc.columns if c.startswith("flag_") and c != "flag_any"]
    aft = pd.DataFrame({
        "criterion": [c.replace("flag_", "") for c in aflags],
        "analytes flagged": [int(aqc[c].sum()) for c in aflags],
        "% of analytes": [_pct(int(aqc[c].sum()), adat.n_analytes) for c in aflags],
        "used for exclusion": ["yes" if c in p["analyte_exclude_on"] else "flag only"
                               for c in aflags]})

    excl = study_sqc[study_sqc.exclude]
    excl_tbl = pd.DataFrame({
        "sample": excl.index.astype(str),
        "plate": excl.get("PlateId", pd.Series("", index=excl.index)).astype(str),
        "reason": excl["exclude_reason"],
        "median log2 RFU": excl["median_log2_rfu"].round(2),
        "r to median": excl["cor_to_median"].round(3)}) if len(excl) else pd.DataFrame(
        columns=["sample", "plate", "reason"])

    prov = pd.DataFrame({"field": [
        "ADAT file", "SHA-256", "AdatId", "Assay version", "Study matrix", "Normalization",
        "SomaLogic process steps", "Wells in file", "Analytes in file", "Grouping column"],
        "value": [
            os.path.basename(p["input_adat"]), adat.sha256[:32] + "...",
            adat.header.get("AdatId", "-"), adat.header.get("AssayVersion", "-"),
            adat.header.get("StudyMatrix", "-"),
            adat.header.get("NormalizationAlgorithm", "-"),
            adat.header.get("ProcessSteps", "-"),
            f"{adat.n_samples} ({', '.join(f'{k}={v}' for k, v in adat.sample_types().items())})",
            f"{adat.n_analytes}", str(p["group_col"] or "none (pooled)")]})

    show = ["group_col", "norm_scale_min", "norm_scale_max", "anml_fraction_min",
            "min_correlation_to_median", "max_intensity_mad_z", "sample_min_detect_rate",
            "keep_analyte_types", "keep_organisms", "cal_qc_ratio_min", "cal_qc_ratio_max",
            "max_replicate_cv", "lod_mad_multiplier", "min_detect_frac", "detect_rule",
            "log_transform", "plate_median_center", "sample_exclude_on", "analyte_exclude_on"]
    ptab = pd.DataFrame({"parameter": show, "value": [str(p[k]) for k in show]})

    gtab = None
    if det_table is not None:
        rows = []
        for c in det_table.columns:
            lvl = c.replace("detect_rate::", "")
            rows.append({
                "group": lvl, "samples": int((groups == lvl).sum()),
                f"analytes detected in >={p['min_detect_frac']:.0%}":
                    int((det_table[c] >= p["min_detect_frac"]).sum()),
                "median detection rate": f"{det_table[c].median():.3f}"})
        gtab = pd.DataFrame(rows)

    otab = pd.DataFrame({
        "property": ["matrix orientation", "proteins (rows)", "samples (columns)",
                     "value scale", "missing values", "value range",
                     "analytes flagged but retained"],
        "value": ["proteins x samples", f"{mat.shape[0]:,}", f"{mat.shape[1]:,}",
                  p["log_transform"], f"{int(mat.isna().values.sum())}",
                  f"{np.nanmin(mat.values):.2f} to {np.nanmax(mat.values):.2f}",
                  f"{int((aqc.loc[keep_analytes, 'flag_any']).sum()):,} "
                  f"(see flag columns in feature metadata)"]})

    lod_cap = ("Buffer wells define the blank distribution; the limit of detection is "
               f"median(buffer) + {p['lod_mad_multiplier']} x MAD(buffer) per analyte."
               if lod is not None else
               "No blank-based limit of detection was estimated for this run.")

    cards = [("study samples in", f"{n_study:,}"), ("samples retained", f"{n_keep:,}"),
             ("analytes in", f"{adat.n_analytes:,}"),
             ("analytes retained", f"{int(keep_analytes.sum()):,}"),
             ("flagged, retained", f"{int(aqc.loc[keep_analytes, 'flag_any'].sum()):,}"),
             ("runtime", "-" if p["deterministic_report"] else f"{elapsed:.0f} s")]

    return {
        "run": {"input_file": os.path.basename(p["input_adat"]),
                "generated": ("(timestamp suppressed: deterministic_report)"
                              if p["deterministic_report"]
                              else _dt.datetime.now().strftime("%Y-%m-%d %H:%M")),
                "version": __version__, "param_hash": param_hash(p),
                "input_sha256": adat.sha256,
                "elapsed_s": (0.0 if p["deterministic_report"] else elapsed)},
        "figures": figs, "summary": {"cards": cards}, "provenance": prov,
        "params_table": ptab, "sample_flag_table": sft, "analyte_flag_table": aft,
        "excluded_samples": excl_tbl, "group_table": gtab, "output_table": otab,
        "lod_caption": lod_cap, "warnings": warnings,
        "output_note": ("Analytes failing a criterion that is set to <em>flag only</em> are kept "
                        "in the matrix and marked in the feature metadata, so downstream "
                        "analyses can exclude or down-weight them without re-running."),
    }
