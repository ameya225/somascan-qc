# QC criteria and methods

This pipeline starts from the ANML-normalized RFU values that SomaLogic delivers and
applies only quality control, filtering and a log transform. It never re-normalizes,
imputes or batch-corrects.

## Why start at ANML

SomaLogic's standard processing chain is hybridization normalization, plate scaling,
calibration, and median signal normalization to a reference (ANML, adaptive normalization
by maximum likelihood). The delivered `...anmlSMP.adat` is the recommended analysis
starting point; re-deriving normalization outside the vendor's calibrator and reference
framework is not reproducible across studies. The pipeline records the ADAT's
`ProcessSteps` header in the report and the manifest so the starting point is explicit.

## Sample-level criteria

| Flag | Default rule | Rationale |
|---|---|---|
| `flag_rowcheck` | `RowCheck` must be `PASS` | SomaLogic's own per-sample acceptance call on normalization. |
| `flag_norm_scale` | every `NormScale_*` in 0.4-2.5 | The median/ANML scale factor per dilution group. A sample needing a scale factor far from 1 differs systematically from the reference - usually volume, handling, or a non-plasma matrix. |
| `flag_hyb_scale` | `HybControlNormScale` in 0.4-2.5 | Hybridization control performance: a readout of the array step rather than the sample. |
| `flag_anml_fraction` | every `ANMLFractionUsed_*` >= 0.3 | Fraction of analytes ANML could match to the reference distribution. A low value means the sample's profile is broadly atypical and its normalization rests on little data. |
| `flag_low_correlation` | Pearson r to the study median log2 profile >= 0.8 | Catches globally aberrant samples that still pass the vendor's checks. |
| `flag_intensity_outlier` | \|MAD-z\| of median log2 RFU <= 5 | Robust outlier test on overall intensity; MAD-z is insensitive to the outliers it is looking for. |
| `flag_low_detection` | fraction of analytes above LoD >= 0.5 | A sample where most analytes sit at the blank level has failed. |
| `flag_sample_note` | flag only | `SampleNotes` is non-empty (e.g. "Low Volume"). Surfaced, never automatic grounds for exclusion. |

Control wells (Buffer, QC, Calibrator) are scored and plotted but never "excluded" -
they are simply not study samples, and they are used to derive limits of detection and
coefficients of variation.

## Analyte-level criteria

| Flag | Default action | Rationale |
|---|---|---|
| `flag_not_protein` | drop | `Type` is not `Protein`: hybridization control elution, spuriomer, non-biotin and non-cleavable controls. These are assay controls, not measurements. |
| `flag_organism` | drop | `Organism` is not in `keep_organisms`. The 11K menu includes non-human reagents; keep them with `--keep-organisms` if your design needs them. |
| `flag_colcheck` | flag only | SomaLogic's `ColCheck` is `FLAG`, meaning the analyte's calibration QC ratio fell outside specification on at least one plate. |
| `flag_cal_qc_ratio` | flag only | Computed independently from the `CalQcRatio_*` columns against the 0.8-1.2 window; `n_plates_calqc_out` records how many plates failed. |
| `flag_high_cv` | flag only | Median intra-plate coefficient of variation across the pooled QC wells exceeds 0.2. Computed within plate and then combined across plates, so a plate offset cannot inflate it. |
| `flag_low_detection` | drop | Detectability, below. |

`ColCheck`/`CalQcRatio` default to flag-only deliberately. Those analytes are measured
but poorly calibrated; whether that matters depends on whether you compare across plates,
so the decision belongs to the analyst and is one parameter away
(`analyte_exclude_on: [..., flag_colcheck]`).

## Limits of detection and the detectability filter

Buffer (blank) wells give the per-analyte noise distribution:

```
LoD = median(buffer RFU) + k * MAD(buffer RFU),    MAD scaled by 1.4826, k = 4.9 by default
```

A protein counts as detected in a sample when its RFU exceeds its LoD. Detection rates
are then computed **within each level of `group_col`**, and the analyte is retained when

- `detect_rule: any` (default) - it reaches `min_detect_frac` in at least one group, so a
  protein present only in one arm of the study survives; or
- `detect_rule: all` - it reaches `min_detect_frac` in every group, giving a matrix with
  no group-dependent missingness at the cost of losing group-specific signal.

`k`, `min_detect_frac` and the rule are parameters. If the ADAT carries fewer than three
buffer wells no LoD can be estimated: the filter is skipped, and both the report and the
manifest record that explicitly rather than silently passing everything.

## What is deliberately not done

- **No imputation.** SomaScan reports a value for every analyte in every sample; values
  below LoD are noise, not missing data, and are kept with their detectability recorded.
- **No batch correction.** `--plate-median-center` exists but is off by default. The
  report shows PCA by plate and residual plate medians so you can decide whether to model
  plate as a covariate instead, which is usually preferable.
- **No re-normalization.** See above.

## Background reading

- Candia J, et al. Assessment of variability in the SomaScan assay. *Sci Rep* 2017;7:14248.
- Candia J, Daya GN, Tanaka T, Ferrucci L, Walker KA. Assessment of variability in the plasma 7k SomaScan proteomics assay. *Sci Rep* 2022;12:17147.
- Joshi A, Mayr M. In aptamers they trust: the caveats of the SOMAscan biomarker discovery platform. *Circulation* 2018;138:2482-2485.
- SomaLogic `SomaDataIO` (CRAN) and the SomaScan assay technical notes, for the ADAT format and the normalization chain.
