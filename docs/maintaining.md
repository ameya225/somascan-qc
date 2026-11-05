# Maintaining and releasing

The repository is published at <https://github.com/ameya225/somascan-qc> with tags
`v1.0.0` and `v1.1.0`.

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request to `main`:

| Job | What it proves |
|---|---|
| tests | pytest on Python 3.10-3.12 (Linux) |
| lint | `ruff check somascan_qc tests` |
| demo | `somascan-qc demo` runs end to end and writes a report, a matrix and a manifest |
| reproducibility | two runs of the same simulated ADAT are compared with `diff`; any byte difference fails the build |
| docker build | the image builds and the demo runs inside it |
| tests (macOS) | **on demand only** - start it from the Actions tab with "Run workflow". Hosted macOS arm64 runners are capacity-constrained and jobs are often never acquired, which ends a run as cancelled; a cancelled job fails the run even with `continue-on-error`, so macOS cannot sit on the push path |

The Docker image and the workflows were never executed in the environment where this
repository was assembled, so the first CI run is their first real test. If the Docker
job fails, the likely causes are the `mambaorg/micromamba` base tag and the
`/opt/conda/bin` paths in the `Dockerfile`.

## Cutting a release

1. Update `CHANGELOG.md`.
2. Bump the version in `somascan_qc/__init__.py`, `pyproject.toml`, `CITATION.cff` and
   `Makefile`.
3. `git tag -a vX.Y.Z -m "vX.Y.Z" && git push --follow-tags`.
4. `gh release create vX.Y.Z --notes-file CHANGELOG.md`, or draft it in the web UI.

## Optional next steps

- **Zenodo DOI**: link the repository on <https://zenodo.org/account/settings/github/>
  before cutting a release; Zenodo mints a DOI for each tagged release. Add the badge to
  the README and the `doi:` field to `CITATION.cff`.
- **Branch protection**: require the `tests` and `docker build` checks on `main`.
- **Topics**: `proteomics`, `somascan`, `somalogic`, `quality-control`, `aptamer`,
  `reproducible-research`.
- **PyPI**: the package builds with `python -m build`; publishing would make
  `pip install somascan-qc` work without a clone.

## What is deliberately not in the repository

No ADAT files, no result matrices, no reports. `.gitignore` excludes `*.adat`, `*.tsv`,
`qc_report.html` and the output directories, so study data cannot be committed by
accident. The test suite and CI run on a simulated ADAT generated at runtime by
`somascan_qc/simulate.py`.
