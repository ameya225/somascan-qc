# Contributing

## Development setup

```bash
git clone https://github.com/ameya225/somascan-qc.git
cd somascan-qc
conda env create -f environment.yml
conda activate somascan-qc
pip install --no-deps -e ".[test]"
make test
```

## Before opening a pull request

```bash
make lint          # ruff
make test          # pytest
make demo          # end-to-end run on a simulated ADAT
```

All three also run in CI on Python 3.10 through 3.12 (Linux). A macOS check is
available on demand from the Actions tab.

## Project conventions

- **No silent data loss.** Every sample and analyte that the pipeline removes must stay
  in `*_sample_qc_all.tsv` / `*_analyte_qc_all.tsv` with a populated `exclude_reason`.
- **New QC criteria** are added as a `flag_<name>` column, are configurable through
  `sample_exclude_on` / `analyte_exclude_on`, default to *flag only* rather than *drop*,
  and get a row in the relevant report table.
- **Every new parameter** goes in `somascan_qc/params.py::DEFAULTS`, in
  `config/params.yaml` with a comment, and in `docs/parameters.md`.
- **Determinism.** No randomness outside `somascan_qc.simulate`. A new feature must not
  break `test_reproducible_outputs`, which asserts identical checksums across two runs.
- **Dependencies** stay limited to numpy, pandas, matplotlib and PyYAML; the report is a
  single self-contained HTML file with no external assets.

## Adding support for another ADAT variant

`somascan_qc/adat.py` discovers the table offset from the file rather than assuming it.
If you meet a file it cannot read, add a generator case to `somascan_qc/simulate.py`
and a parser test before changing the reader.

## Releasing

1. Update `CHANGELOG.md` and bump `__version__` in `somascan_qc/__init__.py`,
   `pyproject.toml` and `CITATION.cff`.
2. `git tag -a vX.Y.Z -m "vX.Y.Z" && git push --tags`.
