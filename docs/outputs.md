# Output files

Every path below is prefixed with `output_prefix` (default `somascan`) and written into
`output_dir`. With `compress_outputs: true` each `.tsv` becomes `.tsv.gz`.

## `*_protein_x_sample_log2.tsv` - the main deliverable

Rows are proteins (`SeqId`), columns are samples. Values are log2 RFU
(`log_transform: log10 | none` changes the scale and the filename suffix). Column names
come from `matrix_column_id` (default: the ADAT's unique sample key, usually
`ExtIdentifier`); `*_sample_metadata.tsv` carries a `matrix_column` column to join on.
There are no missing values.

## `*_protein_x_sample_rfu.tsv`

The same matrix, same filtering, on the untransformed RFU scale.

## `*_feature_metadata.tsv` - one row per retained protein

| Column | Meaning |
|---|---|
| `SeqId` | SOMAmer identifier; the matrix row key |
| `SomaId`, `Target`, `TargetFullName`, `UniProt`, `EntrezGeneSymbol`, `EntrezGeneID` | annotation copied from the ADAT |
| `Organism`, `Type`, `Dilution`, `Units` | analyte class and dilution bin |
| `ColCheck` | SomaLogic's own analyte QC call |
| `flag_not_protein`, `flag_organism`, `flag_colcheck`, `flag_cal_qc_ratio`, `flag_high_cv`, `flag_low_detection` | one boolean per criterion |
| `n_plates_calqc_out` | number of plates whose calibration QC ratio fell outside the window |
| `qc_replicate_cv` | median intra-plate CV across pooled QC wells |
| `lod_rfu` | blank-based limit of detection |
| `detect_rate_overall` | fraction of retained samples above the LoD |
| `detect_rate::<group>` | the same, within each level of `group_col` |
| `n_groups_detected` | how many groups met `min_detect_frac` |
| `flag_any` | any criterion fired |
| `exclude`, `exclude_reason` | always `False` / empty here by construction |

## `*_sample_metadata.tsv` - one row per retained sample

Every ADAT sample-annotation column, plus `group`, `matrix_column`, `median_log2_rfu`,
`cor_to_median`, `detect_rate` and `median_intensity_z`.

## `*_analyte_qc_all.tsv` and `*_sample_qc_all.tsv`

The same tables for **every** analyte and **every** well in the file, including the ones
that were dropped and including control wells, each with `exclude` and a populated
`exclude_reason`. Nothing leaves the pipeline unaccounted for; these are the files to
read when a count looks wrong.

## `qc_report.html`

Self-contained (figures embedded as base64, no external assets, opens offline).
Sections: provenance, parameters, sample QC, signal and limits of detection, analyte QC
and filtering, batch structure, output summary.

## `*_run_manifest.json`

| Key | Contents |
|---|---|
| `input` | absolute path, SHA-256, dimensions, well-type counts, selected ADAT header fields |
| `parameters` | every effective parameter |
| `param_hash` | 12-character digest of the parameters, excluding paths |
| `results` | samples and analytes in/out, exclusions, group sizes |
| `timings_s` | per-step and total wall time. `total_s` is the only complete figure: the report footer is rendered before the run ends, so the number it shows covers everything up to report rendering and is smaller by the report and manifest steps |
| `environment` | Python, platform, numpy and pandas versions |
| `outputs` | path and SHA-256 of every file written |
| `warnings` | anything the run wants you to know (skipped LoD filter, pooled grouping, ...) |

Two runs with the same input and parameters produce identical `sha256` values for every
data output; set `deterministic_report: true` for the HTML to match as well.
