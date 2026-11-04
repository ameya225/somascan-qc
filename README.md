# somascan-qc

[![CI](https://github.com/ameya225/somascan-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/ameya225/somascan-qc/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

Turn a SomaLogic **ADAT** file (ANML-normalized SomaScan 5K / 7K / 11K) into an
**analysis-ready protein x sample matrix** plus a **self-contained HTML QC report**,
in one command.

```
somascan-qc run -i mystudy.adat -o results -g SampleGroup
```

No real data to hand? `somascan-qc demo` simulates an ADAT and runs the whole pipeline
on it, which is also the fastest way to check that an install or container works.

Everything the pipeline decided is written down: every sample and every analyte is
kept in a QC table with one column per criterion, the exact parameters are stored in a
JSON manifest with SHA-256 checksums of every output, and re-running the same input
with the same parameters reproduces byte-identical files.

---

## 1. What you need

* Your `.adat` file, already normalized by SomaLogic (the filename usually ends in
  `...anmlSMP.adat` or `...hybNorm.medNormInt.plateScale.calibrate.anmlQC.anmlSMP.adat`).
  This pipeline **does not re-normalize** — it starts from SomaLogic's delivered
  ANML values, which is the recommended starting point for analysis.
* Either **Docker** (nothing else to install) or **conda**.

## 2. Install

### Option A - Docker (recommended if you are not a bioinformatician)

```bash
cd somascan_pipeline
docker build -t somascan-qc:1.0.0 .
```

Run it (mount the folder that holds your ADAT as `/data`):

```bash
docker run --rm -v "$PWD":/data somascan-qc:1.0.0 \
    run -i /data/mystudy.adat -o /data/results -g SampleGroup
```

### Option B - conda

```bash
cd somascan_pipeline
conda env create -f environment.yml
conda activate somascan-qc
pip install --no-deps -e .
somascan-qc run -i mystudy.adat -o results -g SampleGroup
```

### Option C - one-line wrapper

`run_pipeline.sh` uses the Docker image when it exists and the conda install otherwise:

```bash
./run_pipeline.sh mystudy.adat results SampleGroup
```

Open `results/qc_report.html` in any web browser when it finishes.

## 3. The three commands

| Command | What it does |
|---|---|
| `somascan-qc inspect -i file.adat` | Prints wells, analytes, SomaLogic processing steps and **which columns you can use as `--group-col`**. Run this first. |
| `somascan-qc run -i file.adat -o out` | Runs the pipeline. |
| `somascan-qc params` | Prints every parameter with its default. |

A typical first session:

```bash
somascan-qc inspect -i mystudy.adat          # see what is in the file
somascan-qc run -i mystudy.adat -o results -g SampleGroup
```

## 4. Setting parameters

Anything can be set on the command line or in a config file. A config file is the
reproducible option — copy `config/params.yaml`, edit it, and pass it with `-c`:

```bash
somascan-qc run -c my_params.yaml
```

Command-line flags override the config file. The most used flags:

```
-i, --input-adat        path to the .adat
-o, --output-dir        where results are written
-g, --group-col         ADAT sample column that defines your groups
    --min-detect-frac   fraction of samples a protein must be detected in (default 0.5)
    --detect-rule       any | all  - pass in at least one group, or in every group
    --keep-organisms    e.g. Human Mouse  (give the flag with no values to keep all)
    --no-detectability-filter
    --plate-median-center
    --log-transform     log2 | log10 | none
```

### Grouping

`--group-col` accepts **any populated sample column of the ADAT** - `SampleGroup`,
`TimePoint`, `SiteId`, `PlateId`, or any custom column your study put there.
`somascan-qc inspect` lists the usable ones. If you leave it unset, all samples are
treated as one pooled group and the report says so.

> If your group labels are not in the ADAT, add them as a column to the ADAT's sample
> block (or ask your provider to include them) before running; the pipeline
> deliberately does not join external spreadsheets, so the matrix can never be
> mismatched to a sample sheet.

## 5. What the pipeline does

**Sample (row) QC** — every well is scored; by default a study sample is removed when
any of these fail:

| Criterion | Default | Why |
|---|---|---|
| `RowCheck` | must be PASS | SomaLogic's own per-sample normalization acceptance flag |
| `NormScale_*` | 0.4 - 2.5 | median/ANML scale factor outside the window means the sample needed implausible scaling |
| `HybControlNormScale` | 0.4 - 2.5 | hybridization control performance |
| `ANMLFractionUsed_*` | >= 0.3 | fraction of analytes ANML could match to the reference |
| correlation to study median | >= 0.8 | catches globally aberrant profiles |
| median intensity MAD-z | <= 5 | catches intensity outliers |
| detection rate | >= 0.5 | fraction of analytes above the limit of detection |

**Analyte (column) QC** — every analyte gets one flag column per criterion:

| Flag | Default action | Meaning |
|---|---|---|
| `flag_not_protein` | **drop** | hybridization controls, spuriomers, non-biotin/non-cleavable controls |
| `flag_organism` | **drop** | non-human analytes (keep them with `--keep-organisms`) |
| `flag_low_detection` | **drop** | not detected above the limit of detection in enough samples per group |
| `flag_colcheck` | flag only | SomaLogic `ColCheck` = FLAG (calibration QC ratio out of spec) |
| `flag_cal_qc_ratio` | flag only | calibration QC ratio outside 0.8 - 1.2 on >= 1 plate |
| `flag_high_cv` | flag only | intra-plate CV from the pooled QC wells above 20 % |

**Limit of detection.** Buffer (blank) wells define the noise floor:
`LoD = median(buffer) + 4.9 x MAD(buffer)` per analyte. A protein is "detected" in a
sample when its RFU exceeds its LoD. The detectability filter keeps a protein when its
detection rate reaches `min_detect_frac` in **any** group (`detect_rule: any`, the
default, so group-specific proteins survive) or in **every** group (`detect_rule: all`).
If the ADAT has fewer than 3 buffer wells the filter is skipped and the report says so.

**Flag vs drop.** Which flags actually remove data is a parameter
(`sample_exclude_on`, `analyte_exclude_on`). Everything else is carried into the
feature metadata so you can exclude or down-weight it later without re-running. To be
stricter, add flags to the list:

```yaml
analyte_exclude_on: [flag_not_protein, flag_organism, flag_low_detection,
                     flag_colcheck, flag_high_cv]
```

**Transform.** The matrix is log2(RFU) by default - the conventional scale for
SomaScan analysis. Set `log_transform: none` for raw RFU.

**What the pipeline does *not* do:** it does not re-normalize (ANML is SomaLogic's job),
it does not impute, and it does not batch-correct. Optional per-plate median centering
is available (`--plate-median-center`) but is off by default; the report shows the
residual plate structure so you can decide.

## 6. Outputs

| File | Contents |
|---|---|
| `*_protein_x_sample_log2.tsv` | **the main deliverable** - rows = SeqId (proteins), columns = samples, log2 RFU |
| `*_protein_x_sample_rfu.tsv` | same matrix on the untransformed RFU scale |
| `*_feature_metadata.tsv` | one row per retained protein: SeqId, Target, UniProt, gene symbol, dilution, LoD, per-group detection rates, QC CV and every flag |
| `*_sample_metadata.tsv` | one row per retained sample: all ADAT annotation + group + QC metrics |
| `*_analyte_qc_all.tsv` | every analyte in the file, including the dropped ones, with the reason |
| `*_sample_qc_all.tsv` | every well in the file, including controls, with the reason |
| `qc_report.html` | self-contained report (no internet needed) |
| `*_run_manifest.json` | parameters, input SHA-256, per-step timings, output checksums, software versions |

Reading the matrix:

```r
mat  <- read.delim("results/somascan_protein_x_sample_log2.tsv", row.names = 1, check.names = FALSE)
feat <- read.delim("results/somascan_feature_metadata.tsv", row.names = 1)
```
```python
import pandas as pd
mat  = pd.read_csv("results/somascan_protein_x_sample_log2.tsv", sep="\t", index_col=0)
feat = pd.read_csv("results/somascan_feature_metadata.tsv",      sep="\t", index_col=0)
```

## 7. Reading the QC report

1. **Run and assay provenance** - confirms the file really is ANML-normalized and records its checksum.
2. **Parameters** - exactly what was applied.
3. **Sample-level QC** - scale factors per plate against the acceptance band, ANML fraction, per-sample intensity (excluded samples ringed in red), concordance and detection rate, and a table of every excluded sample with its reason.
4. **Signal and limits of detection** - buffer wells should sit clearly left of the study samples; if they overlap, the assay did not separate signal from noise.
5. **Analyte QC and filtering** - detectability, QC-replicate CV (a healthy SomaScan run has a median CV near 5 %), flag frequencies and the filtering waterfall.
6. **Batch structure** - PCA coloured by plate and by group, and residual plate medians. Plate-driven separation means you should model plate (or enable centering).
7. **Output matrix** - final dimensions and how many retained proteins still carry a flag.

## 8. Reproducibility

* Same input + same parameters = byte-identical `.tsv` outputs (verified by checksum).
  Set `deterministic_report: true` to make the HTML byte-identical too (it otherwise
  carries a generation timestamp).
* `*_run_manifest.json` records the input SHA-256, a 12-character digest of the
  effective parameters, every output checksum, per-step timings and library versions.
* No random number generation is involved anywhere in the pipeline.

## 9. Performance

Measured on a SomaScan 11K ADAT (576 wells x 11,083 analytes, 44 MB, 6 plates) on a
laptop-class 8-core machine (32 GB RAM): **11-12 s end to end**, peak memory 1.3 GB. Roughly 55 % of
that is writing the two 40 MB text matrices - set `compress_outputs: true` to write
`.tsv.gz` instead.

## 10. Tests

```bash
pytest -q          # 13 tests: parser, QC components, end-to-end, determinism
```

The tests build a small synthetic ADAT (with known bad samples, control analytes,
non-human analytes and blank-level analytes) and assert that each is handled correctly,
including that two runs produce identical checksums.

## 11. Troubleshooting

| Message | Fix |
|---|---|
| `not a valid ADAT - missing ^COL_DATA` | The file is not an ADAT, or was opened and re-saved by Excel. Use the original. |
| `group_col 'X' is not a sample column` | Run `somascan-qc inspect` and use one of the listed columns. |
| `every study sample was excluded` | Thresholds are too strict for this run - inspect `*_sample_qc_all.tsv` and relax `sample_exclude_on`. |
| `every analyte was excluded` | Lower `min_detect_frac`, or set `detect_rule: any`, or check the file is ANML-normalized. |
| Report says the LoD filter was skipped | The ADAT has <3 buffer wells. Ask your provider for the version that includes control wells, or set `detectability_filter: false` deliberately. |

## 12. Layout

```
somascan-qc/
├── somascan_qc/
│   ├── adat.py             ADAT reader (header, analyte and sample annotation, RFU)
│   ├── params.py           defaults, validation, parameter digest
│   ├── sampleqc.py         sample-level criteria
│   ├── analyteqc.py        limits of detection, CV, analyte criteria, detectability filter
│   ├── pipeline.py         orchestration, matrix build, outputs, manifest
│   ├── report.py           figures and HTML report
│   ├── simulate.py         synthetic ADAT generator (demo + tests)
│   └── cli.py              command-line interface
├── config/params.yaml      fully commented parameter template
├── docs/                   parameter reference, QC methods, output schemas
├── examples/               ready-made config files
├── tests/                  test suite
├── .github/workflows/ci.yml  tests, demo, reproducibility check, Docker build
├── environment.yml         conda environment (portable)
├── environment.lock.yml    pinned environment (exact tested versions)
├── requirements-lock.txt   pinned pip layer
├── Dockerfile              container image
├── Makefile                make test | lint | demo | docker
└── run_pipeline.sh         one-line wrapper
```

Each module is importable on its own, e.g. to use only the reader:

```python
from somascan_qc.adat import read_adat
adat = read_adat("mystudy.adat")       # adat.rfu, adat.samples, adat.analytes
```

## 13. Documentation

| Document | Contents |
|---|---|
| [`docs/parameters.md`](docs/parameters.md) | every parameter, its default and its meaning |
| [`docs/qc_criteria.md`](docs/qc_criteria.md) | what each QC criterion tests, why, and what the pipeline deliberately does not do |
| [`docs/outputs.md`](docs/outputs.md) | every output file and a column dictionary for each table |
| [`CHANGELOG.md`](CHANGELOG.md) | release history |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | development setup and project conventions |
| [`docs/maintaining.md`](docs/maintaining.md) | CI jobs, release process and optional next steps |

## 14. Citing

If this pipeline contributed to published work, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub renders a "Cite this repository" button from it),
and cite SomaLogic's assay and normalization documentation for the upstream data.

## 15. Licence

MIT - see [`LICENSE`](LICENSE). The pipeline is independent of SomaLogic; "SomaScan",
"SomaLogic" and "SOMAmer" are trademarks of their owner and are used here only to
describe the data format this software reads.
