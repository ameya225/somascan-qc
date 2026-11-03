"""Command-line entry point: ``somascan-qc``."""
from __future__ import annotations

import argparse
import json
import os
import sys

from . import __version__
from .adat import read_adat
from .params import DEFAULTS
from .pipeline import run


def _load_config(path):
    if path is None:
        return {}
    with open(path) as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            sys.exit("PyYAML is required for YAML configs; "
                     "use a .json config or run `pip install pyyaml`")
        return yaml.safe_load(text) or {}
    return json.loads(text)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="somascan-qc",
        description="QC, filter and export a SomaScan ANML-normalized ADAT as a "
                    "protein x sample matrix with an HTML QC report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--version", action="version", version=f"somascan_qc {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the pipeline")
    r.add_argument("-i", "--input-adat", help="path to the .adat file")
    r.add_argument("-o", "--output-dir", help="directory for outputs")
    r.add_argument("-c", "--config", help="YAML/JSON file with parameter overrides")
    r.add_argument("-g", "--group-col", help="ADAT sample column defining groups")
    r.add_argument("--min-detect-frac", type=float,
                   help="min fraction of samples above LoD within a group")
    r.add_argument("--detect-rule", choices=["any", "all"],
                   help="analyte must pass in ANY or ALL groups")
    r.add_argument("--max-replicate-cv", type=float, help="QC-replicate CV flag threshold")
    r.add_argument("--log-transform", choices=["log2", "log10", "none"])
    r.add_argument("--keep-organisms", nargs="*", help="e.g. Human Mouse; empty keeps all")
    r.add_argument("--no-detectability-filter", action="store_true")
    r.add_argument("--plate-median-center", action="store_true")
    r.add_argument("--output-prefix")
    r.add_argument("--quiet", action="store_true")

    i = sub.add_parser("inspect", help="print ADAT structure and usable grouping columns")
    i.add_argument("-i", "--input-adat", required=True)

    sub.add_parser("params", help="print all parameters and their defaults")

    s = sub.add_parser("simulate", help="write a small synthetic ADAT for testing")
    s.add_argument("-o", "--output", default="example.adat")
    s.add_argument("--n-analytes", type=int, default=60)
    s.add_argument("--n-per-plate", type=int, default=24)
    s.add_argument("--seed", type=int, default=0)

    d = sub.add_parser("demo", help="simulate an ADAT and run the pipeline on it")
    d.add_argument("-o", "--output-dir", default="somascan_qc_demo")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.cmd == "simulate":
        from .simulate import write_adat
        truth = write_adat(args.output, n_analytes=args.n_analytes,
                           n_per_plate=args.n_per_plate, seed=args.seed)
        print(f"wrote {args.output}: {truth['n_wells']} wells "
              f"({truth['n_samples']} study samples, {truth['n_buffer']} buffer) "
              f"x {truth['n_analytes']} analytes")
        return 0

    if args.cmd == "demo":
        from .simulate import write_adat
        os.makedirs(args.output_dir, exist_ok=True)
        adat = os.path.join(args.output_dir, "example.adat")
        write_adat(adat)
        man = run({"input_adat": adat, "output_dir": args.output_dir,
                   "group_col": "SampleGroup", "output_prefix": "demo"})
        r = man["results"]
        print(f"demo OK: {r['study_samples']} -> {r['samples_out']} samples, "
              f"{r['analytes_in']} -> {r['analytes_out']} analytes in "
              f"{man['timings_s']['total_s']} s")
        print(f"  report: {man['outputs']['report']['path']}")
        return 0

    if args.cmd == "params":
        print(json.dumps(DEFAULTS, indent=2))
        return 0

    if args.cmd == "inspect":
        a = read_adat(args.input_adat)
        print(f"{a.path}\n  wells    : {a.n_samples} ({a.sample_types()})")
        print(f"  analytes : {a.n_analytes} ({a.analytes['Type'].value_counts().to_dict()})")
        print(f"  steps    : {a.header.get('ProcessSteps','-')}")
        print(f"  norm     : {a.header.get('NormalizationAlgorithm','-')}")
        print("  usable grouping columns (populated, >1 level):")
        for c in a.grouping_candidates():
            lv = a.samples.loc[a.samples.SampleType.eq('Sample'), c].nunique() \
                if 'SampleType' in a.samples else a.samples[c].nunique()
            print(f"    - {c}  ({lv} levels)")
        return 0

    cfg = _load_config(args.config)
    cli = {k: v for k, v in {
        "input_adat": args.input_adat, "output_dir": args.output_dir,
        "group_col": args.group_col, "min_detect_frac": args.min_detect_frac,
        "detect_rule": args.detect_rule, "max_replicate_cv": args.max_replicate_cv,
        "log_transform": args.log_transform, "output_prefix": args.output_prefix,
        "keep_organisms": args.keep_organisms,
        "detectability_filter": False if args.no_detectability_filter else None,
        "plate_median_center": True if args.plate_median_center else None,
    }.items() if v is not None}
    cfg.update(cli)

    try:
        man = run(cfg)
    except (ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")

    if not args.quiet:
        res = man["results"]
        print(f"somascan_qc {man['version']}  ({man['timings_s']['total_s']} s)")
        print(f"  samples : {res['study_samples']} study -> {res['samples_out']} retained "
              f"({res['samples_excluded']} excluded)")
        print(f"  analytes: {res['analytes_in']} -> {res['analytes_out']} retained")
        print(f"  matrix  : {man['outputs']['matrix']['path']}")
        print(f"  report  : {man['outputs']['report']['path']}")
        for w in man["warnings"]:
            print(f"  note    : {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
