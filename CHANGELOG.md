# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [1.1.0] - 2026-11-03

### Added
- `somascan-qc simulate` writes a synthetic ADAT with known defects (bad wells, control
  analytes, non-human analytes, blank-level analytes).
- `somascan-qc demo` simulates an ADAT and runs the full pipeline on it, so an install
  or container can be verified end to end without real study data.
- Repository packaging: MIT licence, citation metadata, contribution guide, changelog,
  GitHub Actions CI (tests on Python 3.9-3.12, end-to-end demo run, Docker build),
  issue and pull-request templates, `Makefile`, pinned dependency lock files and
  extended documentation under `docs/`.

### Changed
- The synthetic-ADAT generator moved from `tests/make_fixture.py` into the installable
  package as `somascan_qc.simulate`, so it ships with the container.
- Lint reformatting under ruff. Data outputs are byte-identical to 1.0.0 on the same
  input (verified by manifest checksums); the HTML report differs only in the version
  string and one CSS line break.

## [1.0.0] - 2026-10-29

### Added
- ADAT reader covering header, analyte annotation, sample annotation and RFU blocks,
  tolerant of the optional spacer column.
- Sample-level QC: `RowCheck`, normalization and hybridization scale-factor windows,
  `ANMLFractionUsed`, correlation to the study median profile, MAD-based intensity
  outliers, per-sample detection rate.
- Analyte-level QC: analyte type, organism, `ColCheck`, calibration QC ratio,
  intra-plate CV from pooled QC wells.
- Blank-based limits of detection from buffer wells and a per-group detectability
  filter with `any`/`all` group rules.
- Configurable flag-versus-drop behaviour for every criterion.
- Outputs: log2 and RFU protein x sample matrices, feature and sample metadata, full
  QC tables for dropped items, a self-contained HTML QC report and a JSON run manifest
  with input and output checksums, per-step timings and library versions.
- Conda environment, Docker image, parameter template and test suite.
