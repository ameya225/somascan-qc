"""Generate a small synthetic ADAT with the same layout as a real SomaLogic file.

Used by the test suite and by ``somascan-qc demo`` so an installation can be verified
end to end without access to real study data. The simulated file contains deliberate
problems - two bad wells, control analytes, non-human analytes and analytes sitting at
the blank level - so that every filter has something to act on.
"""
from __future__ import annotations

import numpy as np

COL_FIELDS = ["SeqId", "SeqIdVersion", "SomaId", "TargetFullName", "Target", "UniProt",
              "EntrezGeneID", "EntrezGeneSymbol", "Organism", "Units", "Type", "Dilution",
              "ColCheck", "CalQcRatio_PLT01_QC1", "CalQcRatio_PLT02_QC1"]
ROW_FIELDS = ["PlateId", "SlideId", "Subarray", "SampleId", "SampleType", "SampleMatrix",
              "SampleNotes", "SampleGroup", "ExtIdentifier", "RowCheck",
              "HybControlNormScale", "NormScale_20", "ANMLFractionUsed_20"]


def write_adat(path, n_analytes=60, n_per_plate=24, spacer=True, seed=0):
    """Write a synthetic ADAT. Returns a dict of the ground truth used by tests."""
    rng = np.random.default_rng(seed)
    plates = ["PLT01", "PLT02"]
    rows = []
    for pl in plates:
        for i in range(n_per_plate):
            if i < 2:
                st, grp = "Buffer", ""
            elif i < 4:
                st, grp = "QC", ""
            else:
                st, grp = "Sample", ("case" if i % 2 else "control")
            rows.append({"PlateId": pl, "SlideId": f"SL{i//4:03d}", "Subarray": str(i % 8 + 1),
                         "SampleId": f"{pl}_{i:02d}", "SampleType": st,
                         "SampleMatrix": "EDTA Plasma", "SampleNotes": "",
                         "SampleGroup": grp, "ExtIdentifier": f"EXID_{pl}_{i:02d}",
                         "RowCheck": "PASS", "HybControlNormScale": "1.0",
                         "NormScale_20": "1.0", "ANMLFractionUsed_20": "0.85"})
    # one sample fails RowCheck, one has an out-of-range normalization scale factor
    bad = [r for r in rows if r["SampleType"] == "Sample"]
    bad[0]["RowCheck"] = "FLAG"
    bad[1]["NormScale_20"] = "3.4"

    seqids = [f"{10000+i}-{i%50}" for i in range(n_analytes)]
    analyte = {f: [] for f in COL_FIELDS}
    for i, sid in enumerate(seqids):
        analyte["SeqId"].append(sid)
        analyte["SeqIdVersion"].append("3")
        analyte["SomaId"].append(f"SL{i:06d}")
        analyte["TargetFullName"].append(f"Protein {i}")
        analyte["Target"].append(f"P{i}")
        analyte["UniProt"].append(f"Q{i:05d}")
        analyte["EntrezGeneID"].append(str(1000 + i))
        analyte["EntrezGeneSymbol"].append(f"GENE{i}")
        analyte["Organism"].append("Mouse" if i % 20 == 0 else "Human")
        analyte["Units"].append("RFU")
        analyte["Type"].append("Spuriomer" if i % 15 == 0 else "Protein")
        analyte["Dilution"].append("20")
        analyte["ColCheck"].append("FLAG" if i % 10 == 0 else "PASS")
        analyte["CalQcRatio_PLT01_QC1"].append("1.30" if i % 10 == 0 else "1.00")
        analyte["CalQcRatio_PLT02_QC1"].append("1.00")

    n = len(rows)
    rfu = np.zeros((n, n_analytes))
    levels = rng.uniform(6.0, 14.0, n_analytes)   # analyte dynamic range, log2 RFU
    for j in range(n_analytes):
        blank_level = (j % 7 == 0)   # every 7th analyte is indistinguishable from buffer
        for i, r in enumerate(rows):
            is_buf = r["SampleType"] == "Buffer"
            if blank_level:
                lvl = levels[j] + (0.5 if is_buf else 0.0)   # samples sit below the blank
            else:
                lvl = levels[j] - (3.0 if is_buf else 0.0)   # real signal, above blank
            rfu[i, j] = 2 ** (lvl + rng.normal(0, 0.15))

    pad = 1 if spacer else 0
    lead = len(ROW_FIELDS) + pad
    with open(path, "w", newline="") as fh:
        def w(s):
            return fh.write(s + "\r\n")
        w("^HEADER")
        for k, v in [("!Version", "1.2"), ("!AssayVersion", "v5.0"),
                     ("!NormalizationAlgorithm", "ANML"), ("!StudyMatrix", "EDTA Plasma"),
                     ("!ProcessSteps", "Raw RFU, hybNorm, plateScale, Calibration, anmlSMP"),
                     ("AdatId", "GID-TEST")]:
            w(f"{k}\t{v}")
        w("^COL_DATA")
        w("!Name\t" + "\t".join(COL_FIELDS))
        w("!Type\t" + "\t".join(["String"] * len(COL_FIELDS)))
        w("^ROW_DATA")
        w("!Name\t" + "\t".join(ROW_FIELDS))
        w("^TABLE_BEGIN")
        for f in COL_FIELDS:
            w("\t" * (lead - 1) + f + "\t" + "\t".join(analyte[f]))
        w("\t".join(ROW_FIELDS) + "\t" * pad + "\t" + "\t".join(seqids))
        for i, r in enumerate(rows):
            vals = [r[f] for f in ROW_FIELDS] + [""] * pad
            w("\t".join(vals) + "\t" + "\t".join(f"{v:.1f}" for v in rfu[i]))
    return {"n_wells": n, "n_analytes": n_analytes,
            "n_samples": sum(r["SampleType"] == "Sample" for r in rows),
            "n_buffer": sum(r["SampleType"] == "Buffer" for r in rows),
            "blank_level_analytes": [seqids[j] for j in range(n_analytes) if j % 7 == 0],
            "mouse": [seqids[j] for j in range(n_analytes) if j % 20 == 0],
            "non_protein": [seqids[j] for j in range(n_analytes) if j % 15 == 0]}
