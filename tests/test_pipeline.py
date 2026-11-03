"""Unit and end-to-end tests. Run with:  pytest -q"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from somascan_qc.adat import read_adat  # noqa: E402
from somascan_qc.analyteqc import limits_of_detection, replicate_cv  # noqa: E402
from somascan_qc.params import merge_params, param_hash  # noqa: E402
from somascan_qc.pipeline import run  # noqa: E402
from somascan_qc.sampleqc import sample_qc  # noqa: E402
from somascan_qc.simulate import write_adat  # noqa: E402


@pytest.fixture(scope="module")
def fixture(tmp_path_factory):
    d = tmp_path_factory.mktemp("adat")
    path = str(d / "test.adat")
    truth = write_adat(path)
    return path, truth


# ---------------------------------------------------------------- parser
def test_parser_dimensions(fixture):
    path, truth = fixture
    a = read_adat(path)
    assert a.n_samples == truth["n_wells"]
    assert a.n_analytes == truth["n_analytes"]
    assert a.sample_types()["Sample"] == truth["n_samples"]
    assert a.header["NormalizationAlgorithm"] == "ANML"
    assert a.rfu.notna().all().all()
    assert list(a.rfu.columns) == list(a.analytes.index)


def test_parser_handles_missing_spacer_column(tmp_path):
    p1, p2 = str(tmp_path / "a.adat"), str(tmp_path / "b.adat")
    write_adat(p1, spacer=True)
    write_adat(p2, spacer=False)
    a, b = read_adat(p1), read_adat(p2)
    assert a.rfu.shape == b.rfu.shape
    assert np.allclose(a.rfu.values, b.rfu.values)
    assert list(a.analytes.Target) == list(b.analytes.Target)


def test_parser_rejects_non_adat(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("not an adat\n")
    with pytest.raises(ValueError):
        read_adat(str(p))


def test_grouping_candidates(fixture):
    a = read_adat(fixture[0])
    assert "SampleGroup" in a.grouping_candidates()
    assert "ExtIdentifier" not in a.grouping_candidates()


# ------------------------------------------------------------ components
def test_lod_formula():
    rng = np.random.default_rng(1)
    buf = pd.DataFrame(rng.normal(100, 10, (20, 5)), columns=list("abcde"))
    lod = limits_of_detection(buf, 4.9)
    med = buf.median()
    mad = (buf - med).abs().median() * 1.4826
    assert np.allclose(lod.values, (med + 4.9 * mad).values)


def test_replicate_cv_is_intra_plate():
    rfu = pd.DataFrame({"a": [10, 11, 10.5, 100, 101, 100.5]})
    plate = pd.Series(["P1"] * 3 + ["P2"] * 3, index=rfu.index)
    cv_within = replicate_cv(rfu, plate)["a"]
    cv_pooled = replicate_cv(rfu, None)["a"]
    assert cv_within < 0.1 < cv_pooled      # plate offset must not inflate the CV


def test_sample_qc_flags_bad_wells(fixture):
    a = read_adat(fixture[0])
    p = merge_params({})
    qc = sample_qc(a.samples, a.rfu, p)
    assert qc.flag_rowcheck.sum() == 1
    assert qc.flag_norm_scale.sum() == 1
    assert qc.loc[qc.SampleType.eq("Sample"), "exclude"].sum() == 2


def test_param_validation():
    with pytest.raises(ValueError):
        merge_params({"detect_rule": "sometimes"})
    with pytest.raises(ValueError):
        merge_params({"norm_scale_min": 3.0, "norm_scale_max": 1.0})
    with pytest.raises(ValueError):
        merge_params({"not_a_parameter": 1})
    assert param_hash(merge_params({})) == param_hash(merge_params({"output_dir": "x"}))
    assert param_hash(merge_params({})) != param_hash(merge_params({"min_detect_frac": 0.9}))


# -------------------------------------------------------------- pipeline
def test_end_to_end(fixture, tmp_path):
    path, truth = fixture
    out = str(tmp_path / "out")
    man = run({"input_adat": path, "output_dir": out, "group_col": "SampleGroup"})
    for key in ("matrix", "feature_metadata", "sample_metadata", "report",
                "analyte_qc_all", "sample_qc_all"):
        assert os.path.exists(man["outputs"][key]["path"]), key
    mat = pd.read_csv(man["outputs"]["matrix"]["path"], sep="\t", index_col=0)
    feat = pd.read_csv(man["outputs"]["feature_metadata"]["path"], sep="\t", index_col=0)
    smeta = pd.read_csv(man["outputs"]["sample_metadata"]["path"], sep="\t", index_col=0)

    assert mat.shape == (len(feat), len(smeta))                 # proteins x samples
    assert man["results"]["samples_out"] == truth["n_samples"] - 2
    assert set(mat.index) == set(feat.index)
    # control analytes, mouse analytes and blank-level analytes must be gone
    assert not set(truth["non_protein"]) & set(mat.index)
    assert not set(truth["mouse"]) & set(mat.index)
    assert not set(truth["blank_level_analytes"]) & set(mat.index)
    # ColCheck FLAG is flag-only by default: retained, but marked
    assert feat.flag_colcheck.sum() > 0
    assert mat.notna().all().all()
    assert 0 < mat.values.mean() < 20                           # log2 scale
    with open(man["outputs"]["report"]["path"]) as fh:
        html = fh.read()
    assert "SomaScan QC report" in html and "data:image/png;base64," in html
    assert set(smeta["group"]) == {"case", "control"}


def test_detect_rule_all_is_stricter(fixture, tmp_path):
    path, _ = fixture
    a = run({"input_adat": path, "output_dir": str(tmp_path / "any"),
             "group_col": "SampleGroup", "detect_rule": "any", "min_detect_frac": 0.9})
    b = run({"input_adat": path, "output_dir": str(tmp_path / "all"),
             "group_col": "SampleGroup", "detect_rule": "all", "min_detect_frac": 0.9})
    assert b["results"]["analytes_out"] <= a["results"]["analytes_out"]


def test_flag_only_vs_exclude(fixture, tmp_path):
    path, _ = fixture
    base = run({"input_adat": path, "output_dir": str(tmp_path / "b")})
    strict = run({"input_adat": path, "output_dir": str(tmp_path / "s"),
                  "analyte_exclude_on": ["flag_not_protein", "flag_organism",
                                         "flag_low_detection", "flag_colcheck"]})
    assert strict["results"]["analytes_out"] < base["results"]["analytes_out"]


def test_bad_group_col_raises(fixture, tmp_path):
    with pytest.raises(ValueError, match="not a sample column"):
        run({"input_adat": fixture[0], "output_dir": str(tmp_path / "x"),
             "group_col": "NoSuchColumn"})


def test_reproducible_outputs(fixture, tmp_path):
    path, _ = fixture
    cfg = {"input_adat": path, "group_col": "SampleGroup"}
    m1 = run({**cfg, "output_dir": str(tmp_path / "r1")})
    m2 = run({**cfg, "output_dir": str(tmp_path / "r2")})
    data_keys = [k for k in m1["outputs"] if k != "report"]
    assert all(m1["outputs"][k]["sha256"] == m2["outputs"][k]["sha256"] for k in data_keys)
    assert m1["param_hash"] == m2["param_hash"]
    # the HTML differs only in its timestamp line
    with open(m1["outputs"]["report"]["path"]) as fh:
        h1 = fh.read()
    with open(m2["outputs"]["report"]["path"]) as fh:
        h2 = fh.read()
    assert len(h1) == len(h2)
