"""Reader for SomaLogic ADAT files (SomaScan v4.x/v5.x, 5K/7K/11K menus).

The ADAT is a tab-delimited text file with four sections:

    ^HEADER        key/value run metadata (assay version, normalization steps, ...)
    ^COL_DATA      names/types of the per-ANALYTE annotation fields
    ^ROW_DATA      names of the per-SAMPLE annotation fields
    ^TABLE_BEGIN   one line per COL_DATA field (analyte annotation, transposed),
                   then a header line of the ROW_DATA field names,
                   then one line per sample: sample annotation + RFU values.

Analyte-annotation lines carry the field name in a spacer column placed after the
sample-annotation block; the RFU block starts one column further right. The offset
is discovered from the file rather than assumed, so files with or without the
spacer column parse identically.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import pandas as pd


def _sha256(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class Adat:
    """Parsed ADAT: header metadata, analyte annotation, sample annotation, RFU."""

    header: Dict[str, str]
    analytes: pd.DataFrame  # index = SeqId, columns = COL_DATA fields
    samples: pd.DataFrame   # index = sample key, columns = ROW_DATA fields
    rfu: pd.DataFrame       # index = sample key, columns = SeqId
    path: str
    sha256: str

    @property
    def n_samples(self) -> int:
        return self.rfu.shape[0]

    @property
    def n_analytes(self) -> int:
        return self.rfu.shape[1]

    def sample_types(self) -> Dict[str, int]:
        if "SampleType" not in self.samples:
            return {}
        return self.samples["SampleType"].value_counts().to_dict()

    def subset_samples(self, mask) -> "Adat":
        mask = np.asarray(mask, dtype=bool)
        return Adat(self.header, self.analytes, self.samples.loc[mask].copy(),
                    self.rfu.loc[mask].copy(), self.path, self.sha256)

    def grouping_candidates(self, max_levels: int = 50) -> List[str]:
        """Sample-annotation columns usable as a grouping variable.

        A column qualifies when it is fully populated across study samples,
        has more than one level, and does not look like a per-sample identifier.
        """
        skip = {"sample_key", "RowCheck", "ColCheck", "Barcode", "Barcode2d",
                "SlideId", "Subarray", "PlatePosition", "ScannerID", "ExtIdentifier",
                "SampleId", "SampleType", "PercentDilution", "SampleNotes"}
        out = []
        st = self.samples.get("SampleType")
        sub = self.samples if st is None else self.samples.loc[st == "Sample"]
        if len(sub) == 0:
            sub = self.samples
        for col in sub.columns:
            if col in skip or col.startswith(("NormScale", "ANMLFraction", "HybControl")):
                continue
            vals = sub[col].astype(str).str.strip()
            if (vals == "").any():
                continue
            n = vals.nunique()
            if 1 < n <= min(max_levels, max(2, len(sub) // 2)):
                out.append(col)
        return out


def read_adat(path: str) -> Adat:
    """Parse an ADAT file into an :class:`Adat`."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]

    sections = {}
    for i, ln in enumerate(lines):
        if ln.startswith("^"):
            sections[ln.split("\t")[0].strip()] = i
    for req in ("^COL_DATA", "^ROW_DATA", "^TABLE_BEGIN"):
        if req not in sections:
            raise ValueError(f"{path}: not a valid ADAT - missing {req} section")

    # ---- header -------------------------------------------------------
    header: Dict[str, str] = {}
    for ln in lines[sections.get("^HEADER", -1) + 1: sections["^COL_DATA"]]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        header[parts[0].lstrip("!").strip()] = parts[1].strip() if len(parts) > 1 else ""

    def _named_row(start: int, end: int, key: str) -> List[str]:
        for ln in lines[start:end]:
            f = ln.split("\t")
            if f[0].strip().lstrip("!").lower() == key:
                return [x.strip() for x in f[1:]]
        raise ValueError(f"{path}: could not find !{key} row")

    col_fields = _named_row(sections["^COL_DATA"] + 1, sections["^ROW_DATA"], "name")
    row_fields = _named_row(sections["^ROW_DATA"] + 1, sections["^TABLE_BEGIN"], "name")

    tb = sections["^TABLE_BEGIN"]
    meta_lines = lines[tb + 1: tb + 1 + len(col_fields)]

    # Discover where the RFU block starts: the analyte-annotation lines carry the
    # COL_DATA field name in a spacer column just left of the first RFU column.
    first = meta_lines[0].split("\t")
    try:
        name_idx = first.index(col_fields[0])
    except ValueError as exc:  # pragma: no cover - malformed file
        raise ValueError(f"{path}: cannot locate analyte annotation offset") from exc
    off = name_idx + 1

    analytes = pd.DataFrame(
        {col_fields[k]: meta_lines[k].split("\t")[off:] for k in range(len(col_fields))}
    )
    if "SeqId" not in analytes.columns:
        raise ValueError(f"{path}: ADAT has no SeqId analyte annotation")
    analytes.index = analytes["SeqId"].astype(str)
    analytes.index.name = "SeqId"

    body = [ln for ln in lines[tb + 1 + len(col_fields):]
            if ln.strip() and not ln.startswith("^")]
    hdr = body[0].split("\t")
    if [x.strip() for x in hdr[:len(row_fields)]] != row_fields:
        raise ValueError(f"{path}: sample annotation header does not match ^ROW_DATA")
    rows = [ln.split("\t") for ln in body[1:]]
    if not rows:
        raise ValueError(f"{path}: ADAT contains no samples")

    samples = pd.DataFrame([r[: len(row_fields)] for r in rows], columns=row_fields)
    rfu = pd.DataFrame([r[off:] for r in rows]).apply(pd.to_numeric, errors="coerce")
    if rfu.shape[1] != analytes.shape[0]:
        raise ValueError(
            f"{path}: {rfu.shape[1]} RFU columns but {analytes.shape[0]} analytes")
    rfu.columns = analytes.index

    key = None
    for cand in ("ExtIdentifier", "SampleId", "Barcode"):
        if cand in samples.columns and samples[cand].astype(str).str.strip().ne("").all() \
                and samples[cand].is_unique:
            key = cand
            break
    idx = samples[key].astype(str) if key else pd.Series(
        [f"row{i+1}" for i in range(len(samples))])
    samples.insert(0, "sample_key", idx.values)
    samples.index = pd.Index(idx.values, name="sample_key")
    rfu.index = samples.index

    return Adat(header=header, analytes=analytes, samples=samples, rfu=rfu,
                path=str(path), sha256=_sha256(path))
