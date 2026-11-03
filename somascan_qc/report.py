"""Self-contained HTML QC report (figures embedded as base64 PNG)."""
from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = ["#2b6cb0", "#dd6b20", "#2f855a", "#9f7aea", "#c53030",
           "#319795", "#b7791f", "#4a5568", "#d53f8c", "#4c51bf"]
FAIL = "#c53030"


def _style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 130, "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "figure.autolayout": False,
    })


def _png(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _band(ax, lo, hi, label=None):
    ax.axhspan(lo, hi, color="#38a169", alpha=0.08, zorder=0)
    for y in (lo, hi):
        ax.axhline(y, color="#38a169", lw=0.8, ls="--", zorder=1)
    if label:
        ax.text(0.99, hi, label, transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=7, color="#2f855a")


# --------------------------------------------------------------------------- figures
def fig_norm_scale(sqc: pd.DataFrame, p: dict) -> str:
    cols = [c for c in sqc.columns if c.startswith("NormScale") or c.startswith("HybControl")]
    if not cols:
        return ""
    study = sqc[sqc.SampleType.eq("Sample")]
    plates = sorted(study["PlateId"].unique()) if "PlateId" in study else ["all"]
    fig, axes = plt.subplots(1, len(cols), figsize=(3.1 * len(cols), 3.3), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, c in zip(axes, cols):
        data = [study.loc[study.PlateId.eq(pl), c].dropna().values if "PlateId" in study
                else study[c].dropna().values for pl in plates]
        bp = ax.boxplot(
            data, patch_artist=True, widths=0.6,
            flierprops={"marker": "o", "markersize": 2,
                        "markerfacecolor": "#718096", "markeredgecolor": "none"})
        for i, box in enumerate(bp["boxes"]):
            box.set(facecolor=PALETTE[i % len(PALETTE)], alpha=0.45, edgecolor="#2d3748")
        for med in bp["medians"]:
            med.set(color="#1a202c", lw=1.2)
        _band(ax, p["norm_scale_min"], p["norm_scale_max"], "acceptance")
        ax.set_title(c, fontsize=9)
        ax.set_xticklabels([str(x)[-5:] for x in plates], rotation=45, ha="right", fontsize=7)
    axes[0].set_ylabel("normalization scale factor")
    fig.suptitle("Normalization scale factors by plate (study samples)", y=1.02, fontsize=11)
    return _png(fig)


def fig_anml(sqc: pd.DataFrame, p: dict) -> str:
    cols = [c for c in sqc.columns if c.startswith("ANMLFractionUsed")]
    if not cols:
        return ""
    study = sqc[sqc.SampleType.eq("Sample")]
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    bins = np.linspace(0, 1, 41)
    for i, c in enumerate(cols):
        ax.hist(study[c].dropna(), bins=bins, histtype="step", lw=1.6,
                color=PALETTE[i % len(PALETTE)], label=c.replace("ANMLFractionUsed_", "dil "))
    ax.axvline(p["anml_fraction_min"], color=FAIL, ls="--", lw=1.1)
    ax.text(p["anml_fraction_min"], ax.get_ylim()[1] * 0.95, f"  min {p['anml_fraction_min']}",
            color=FAIL, fontsize=8, va="top")
    ax.set_xlabel("ANML fraction used")
    ax.set_ylabel("samples")
    ax.legend(fontsize=8)
    ax.set_title("Fraction of analytes used by ANML normalization")
    return _png(fig)


def fig_sample_intensity(sqc: pd.DataFrame) -> str:
    study = sqc[sqc.SampleType.eq("Sample")].copy()
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    if "PlateId" in study:
        plates = sorted(study.PlateId.unique())
        x0 = 0
        for i, pl in enumerate(plates):
            sub = study[study.PlateId.eq(pl)]
            x = np.arange(x0, x0 + len(sub))
            ax.scatter(x, sub.median_log2_rfu, s=7, color=PALETTE[i % len(PALETTE)],
                       label=str(pl), alpha=0.8, edgecolors="none")
            bad = sub[sub.exclude]
            if len(bad):
                ax.scatter(x[sub.exclude.values], bad.median_log2_rfu, s=42, facecolors="none",
                           edgecolors=FAIL, lw=1.2, zorder=3)
            x0 += len(sub) + 4
        ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  title="plate", title_fontsize=7.5)
    else:
        ax.scatter(np.arange(len(study)), study.median_log2_rfu, s=7, color=PALETTE[0])
    ax.set_xlabel("sample (ordered within plate)")
    ax.set_ylabel("median log2 RFU")
    ax.set_title("Per-sample median intensity; red rings = excluded samples")
    return _png(fig)


def fig_sample_hists(sqc: pd.DataFrame, p: dict) -> str:
    study = sqc[sqc.SampleType.eq("Sample")]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
    axes[0].hist(study.cor_to_median.dropna(), bins=40, color=PALETTE[0], alpha=0.8)
    axes[0].axvline(p["min_correlation_to_median"], color=FAIL, ls="--", lw=1.1)
    axes[0].set_xlabel("Pearson r to study median profile")
    axes[0].set_ylabel("samples")
    axes[0].set_title("Global sample concordance")
    if study.detect_rate.notna().any():
        axes[1].hist(study.detect_rate.dropna(), bins=40, color=PALETTE[2], alpha=0.8)
        axes[1].axvline(p["sample_min_detect_rate"], color=FAIL, ls="--", lw=1.1)
        axes[1].set_xlabel("fraction of analytes > LoD")
        axes[1].set_title("Per-sample detection rate")
    else:
        axes[1].axis("off")
    axes[1].set_ylabel("samples")
    return _png(fig)


def fig_intensity_by_type(rfu: pd.DataFrame, stype: pd.Series) -> str:
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    for i, t in enumerate([t for t in ["Buffer", "Calibrator", "QC", "Sample"]
                           if (stype == t).any()]):
        vals = np.log2(rfu.loc[stype.eq(t).values].clip(lower=1e-6).values.ravel())
        ax.hist(vals, bins=120, density=True, histtype="step", lw=1.5,
                color=PALETTE[i % len(PALETTE)], label=f"{t} (n={(stype==t).sum()})")
    ax.set_xlabel("log2 RFU")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    ax.set_title("Signal distribution by well type")
    return _png(fig)


def fig_pca(log2_mat: pd.DataFrame, plate: pd.Series, group: pd.Series) -> str:
    if log2_mat.shape[0] < 3 or log2_mat.shape[1] < 2:
        return ""
    X = log2_mat.values.astype(float)
    X = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=1)
    sd[sd == 0] = 1
    X = X / sd
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pc = U[:, :2] * S[:2]
    var = (S ** 2) / (S ** 2).sum() * 100
    panels = [("plate", plate)]
    if group is not None and group.nunique() > 1 and not group.equals(plate):
        panels.append(("group", group))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.6 * len(panels), 3.8))
    axes = np.atleast_1d(axes)
    for ax, (name, lab) in zip(axes, panels):
        lab = lab.reindex(log2_mat.index).astype(str)
        for i, lv in enumerate(sorted(lab.unique())):
            m = (lab == lv).values
            ax.scatter(pc[m, 0], pc[m, 1], s=11, alpha=0.8, edgecolors="none",
                       color=PALETTE[i % len(PALETTE)], label=lv)
        ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
        ax.set_title(f"PCA coloured by {name}")
        if lab.nunique() <= 10:
            ax.legend(fontsize=7, markerscale=1.4)
    fig.suptitle("PCA of retained analytes (log2, z-scored)", y=1.03, fontsize=11)
    return _png(fig)


def fig_analyte_panels(aqc: pd.DataFrame, p: dict) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.1))
    prot = aqc[~aqc.flag_not_protein]
    if prot.detect_rate_overall.notna().any():
        axes[0].hist(prot.detect_rate_overall.dropna(), bins=50, color=PALETTE[2], alpha=0.85)
        axes[0].axvline(p["min_detect_frac"], color=FAIL, ls="--", lw=1.1)
        axes[0].set_yscale("log")
    axes[0].set_xlabel("fraction of samples > LoD")
    axes[0].set_ylabel("analytes (log scale)")
    axes[0].set_title("Detectability")
    if prot.qc_replicate_cv.notna().any():
        cv = prot.qc_replicate_cv.dropna()
        axes[1].hist(cv, bins=np.linspace(0, max(0.6, float(cv.quantile(0.995))), 60),
                     color=PALETTE[1], alpha=0.85)
        axes[1].axvline(p["max_replicate_cv"], color=FAIL, ls="--", lw=1.1)
        axes[1].set_title(f"QC-replicate CV (median {cv.median():.3f})")
    else:
        axes[1].set_title("QC-replicate CV (no QC wells)")
    axes[1].set_xlabel("intra-plate CV")
    axes[1].set_ylabel("analytes")
    flags = [c for c in aqc.columns if c.startswith("flag_") and c != "flag_any"]
    counts = aqc[flags].sum().sort_values()
    axes[2].barh([f.replace("flag_", "") for f in counts.index], counts.values,
                 color=PALETTE[0], alpha=0.85)
    for y, v in enumerate(counts.values):
        axes[2].text(v, y, f" {int(v)}", va="center", fontsize=7.5)
    axes[2].set_xlabel("analytes flagged")
    axes[2].set_title("Flag frequency")
    axes[2].margins(x=0.18)
    return _png(fig)


def fig_waterfall(steps: list[tuple]) -> str:
    labels = [s[0] for s in steps][::-1]
    vals = [s[1] for s in steps][::-1]
    fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(labels) + 1.4))
    ax.barh(labels, vals, color=PALETTE[0], alpha=0.85)
    for y, v in enumerate(vals):
        ax.text(v, y, f" {v:,}", va="center", fontsize=8)
    ax.set_xlabel("analytes retained")
    ax.set_title("Analyte filtering waterfall")
    ax.margins(x=0.16)
    return _png(fig)


def fig_plate_effect(log2_mat: pd.DataFrame, plate: pd.Series) -> str:
    plate = plate.reindex(log2_mat.index).astype(str)
    if plate.nunique() < 2:
        return ""
    fig, ax = plt.subplots(figsize=(5.8, 3.2))
    med = log2_mat.median(axis=1)
    plates = sorted(plate.unique())
    data = [med[(plate == pl).values].values for pl in plates]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6)
    for i, box in enumerate(bp["boxes"]):
        box.set(facecolor=PALETTE[i % len(PALETTE)], alpha=0.45, edgecolor="#2d3748")
    for m in bp["medians"]:
        m.set(color="#1a202c", lw=1.2)
    ax.set_xticklabels(plates, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("median log2 RFU")
    ax.set_title("Residual plate effect after filtering")
    return _png(fig)


# --------------------------------------------------------------------------- html
_CSS = """
body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
margin:0;color:#1a202c;background:#f7fafc}
.wrap{max-width:1060px;margin:0 auto;padding:28px 32px 64px}
h1{font-size:25px;margin:0 0 4px} h2{font-size:18px;margin:34px 0 10px;
border-bottom:2px solid #e2e8f0;padding-bottom:5px} h3{font-size:14px;margin:18px 0 6px}
.sub{color:#4a5568;font-size:13px;margin-bottom:18px}
table{border-collapse:collapse;font-size:12.5px;margin:10px 0;background:#fff}
th,td{border:1px solid #e2e8f0;padding:5px 10px;text-align:left}
th{background:#edf2f7;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
img{max-width:100%;margin:10px 0;background:#fff;border:1px solid #e2e8f0;border-radius:4px}
.cards{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0}
.card{background:#fff;border:1px solid #e2e8f0;border-radius:6px;padding:12px 16px;min-width:150px}
.card .v{font-size:22px;font-weight:600} .card .l{font-size:11.5px;color:#4a5568;margin-top:2px}
.pass{color:#2f855a} .warn{color:#b7791f} .fail{color:#c53030}
code{background:#edf2f7;padding:1px 5px;border-radius:3px;font-size:12px}
.note{background:#fffaf0;border-left:3px solid #dd6b20;padding:9px 13px;
font-size:12.5px;margin:12px 0}
footer{margin-top:40px;font-size:11.5px;color:#718096;border-top:1px solid #e2e8f0;padding-top:12px}
"""


def _tbl(df: pd.DataFrame, num_cols=()) -> str:
    th = "".join(f"<th>{c}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        tds = "".join(
            f'<td class="num">{r[c]}</td>' if c in num_cols else f"<td>{r[c]}</td>"
            for c in df.columns)
        rows.append(f"<tr>{tds}</tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def _img(b64: str, cap: str = "") -> str:
    if not b64:
        return ""
    c = f'<div class="sub">{cap}</div>' if cap else ""
    return f'<img src="data:image/png;base64,{b64}"/>{c}'


def build_report(ctx: dict) -> str:
    """Assemble the HTML QC report from a context dict built by the pipeline."""
    f = ctx["figures"]
    s = ctx["summary"]
    cards = "".join(
        f'<div class="card"><div class="v">{v}</div><div class="l">{l}</div></div>'
        for l, v in s["cards"])
    parts = [f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>SomaScan QC report - {ctx['run']['input_file']}</title><style>{_CSS}</style></head><body>
<div class="wrap"><h1>SomaScan QC report</h1>
<div class="sub">{ctx['run']['input_file']} &middot; generated {ctx['run']['generated']}
 &middot; somascan_qc v{ctx['run']['version']}</div>
<div class="cards">{cards}</div>"""]

    parts.append("<h2>1. Run and assay provenance</h2>")
    parts.append(_tbl(ctx["provenance"], num_cols=()))
    if ctx["warnings"]:
        parts.append("".join(f'<div class="note">{w}</div>' for w in ctx["warnings"]))

    parts.append("<h2>2. Parameters</h2>")
    parts.append(_tbl(ctx["params_table"]))

    parts.append("<h2>3. Sample-level QC</h2>")
    parts.append(_tbl(ctx["sample_flag_table"], num_cols=("samples flagged", "% of study samples")))
    parts.append(_img(f.get("norm_scale"),
                      "Scale factors outside the acceptance window indicate a sample whose "
                      "normalization is not supported by the data."))
    parts.append(_img(f.get("anml")))
    parts.append(_img(f.get("sample_intensity")))
    parts.append(_img(f.get("sample_hists")))
    if ctx["excluded_samples"] is not None and len(ctx["excluded_samples"]):
        parts.append("<h3>Excluded samples</h3>")
        parts.append(_tbl(ctx["excluded_samples"]))

    parts.append("<h2>4. Signal, controls and limits of detection</h2>")
    parts.append(_img(f.get("intensity_by_type"), ctx["lod_caption"]))

    parts.append("<h2>5. Analyte-level QC and filtering</h2>")
    parts.append(_tbl(ctx["analyte_flag_table"], num_cols=("analytes flagged", "% of analytes")))
    parts.append(_img(f.get("analyte_panels")))
    parts.append(_img(f.get("waterfall")))
    if ctx["group_table"] is not None:
        parts.append("<h3>Detectability by group</h3>")
        parts.append(_tbl(ctx["group_table"]))

    parts.append("<h2>6. Batch structure of the filtered matrix</h2>")
    parts.append(_img(f.get("pca")))
    parts.append(_img(f.get("plate_effect")))

    parts.append("<h2>7. Output matrix</h2>")
    parts.append(_tbl(ctx["output_table"]))
    parts.append(f'<div class="note">{ctx["output_note"]}</div>')

    parts.append(f"""<footer>Pipeline: somascan_qc v{ctx['run']['version']} &middot;
 config digest <code>{ctx['run']['param_hash']}</code> &middot;
 input SHA-256 <code>{ctx['run']['input_sha256'][:16]}...</code> &middot;
 wall time {ctx['run']['elapsed_s']:.1f} s</footer></div></body></html>""")
    return "\n".join(parts)
