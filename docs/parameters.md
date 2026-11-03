# Parameter reference

Every parameter can be set in a YAML/JSON config (`-c`) or, for the common ones, as a
command-line flag. Flags override the config file. `somascan-qc params` prints the
defaults as JSON; `config/params.yaml` is a commented template.

## input/output

| Parameter | Default | Meaning |
|---|---|---|
| `input_adat` | `null` | Path to the ANML-normalized `.adat` file. Required. |
| `output_dir` | `somascan_qc_out` | Directory for all outputs; created if absent. |
| `sample_types_to_keep` | `['Sample']` | ADAT `SampleType` values that enter the matrix. Control wells are still used for LoD and CV. |
| `output_prefix` | `somascan` | Prefix for every output filename. |

## grouping

| Parameter | Default | Meaning |
|---|---|---|
| `group_col` | `null` | Any populated ADAT sample column defining groups for the detectability filter. `null` pools all samples. `somascan-qc inspect` lists the usable columns. |

## sample QC

| Parameter | Default | Meaning |
|---|---|---|
| `norm_scale_min` | `0.4` | Lower bound of the acceptance window for every `NormScale_*` column. |
| `norm_scale_max` | `2.5` | Upper bound of the same window. |
| `hyb_scale_min` | `0.4` | Lower bound for `HybControlNormScale`. |
| `hyb_scale_max` | `2.5` | Upper bound for `HybControlNormScale`. |
| `anml_fraction_min` | `0.3` | Minimum `ANMLFractionUsed_*`: fraction of analytes ANML matched to the reference. |
| `min_correlation_to_median` | `0.8` | Minimum Pearson r between a sample's log2 profile and the study median profile. |
| `max_intensity_mad_z` | `5.0` | Maximum absolute MAD-z of the per-sample median log2 RFU. |
| `sample_min_detect_rate` | `0.5` | Minimum fraction of analytes above their LoD in a sample. |
| `sample_exclude_on` | `['flag_rowcheck', 'flag_norm_scale', 'flag_hyb_scale', 'flag_anml_fraction', 'flag_low_correlation', 'flag_intensity_outlier', 'flag_low_detection']` | Which `flag_*` columns remove a sample. Flags not listed are reported only. |

## analyte QC

| Parameter | Default | Meaning |
|---|---|---|
| `keep_analyte_types` | `['Protein']` | ADAT analyte `Type` values to keep; everything else (hybridization controls, spuriomers, non-biotin, non-cleavable) is flagged. |
| `keep_organisms` | `['Human']` | Organisms to keep. `[]` disables the filter and keeps non-human analytes. |
| `cal_qc_ratio_min` | `0.8` | Lower bound of the calibration QC ratio acceptance window. |
| `cal_qc_ratio_max` | `1.2` | Upper bound of the same window. |
| `max_replicate_cv` | `0.2` | Intra-plate CV from pooled QC wells above which an analyte is flagged. |
| `detectability_filter` | `True` | Enable the blank-based detectability filter. Skipped automatically with fewer than 3 buffer wells. |
| `lod_mad_multiplier` | `4.9` | `k` in `LoD = median(buffer) + k * MAD(buffer)`. |
| `min_detect_frac` | `0.5` | Fraction of samples within a group in which an analyte must exceed its LoD. |
| `detect_rule` | `any` | `any` keeps an analyte passing in at least one group; `all` requires every group. |
| `analyte_exclude_on` | `['flag_not_protein', 'flag_organism', 'flag_low_detection']` | Which `flag_*` columns remove an analyte. Flags not listed are reported only. |

## transform/output

| Parameter | Default | Meaning |
|---|---|---|
| `log_transform` | `log2` | `log2` (default), `log10` or `none`. |
| `write_rfu_matrix` | `True` | Also write the untransformed RFU matrix. |
| `compress_outputs` | `False` | Write `.tsv.gz` instead of `.tsv`. |
| `matrix_column_id` | `sample_key` | ADAT sample column used for the matrix column names; must be unique across retained samples. |
| `plate_median_center` | `False` | Shift each analyte's per-plate median onto its overall median. Off by default. |
| `report_filename` | `qc_report.html` | Filename of the HTML report inside `output_dir`. |
| `deterministic_report` | `False` | Omit the timestamp and runtime so the HTML is byte-identical between runs. |

## Validation

Unknown keys, an inverted scale-factor window, a `detect_rule` other than `any`/`all`,
a `log_transform` other than `log2`/`log10`/`none`, a non-list where a list is
required, or a `group_col` / `matrix_column_id` that is not a sample column of the
ADAT all fail immediately with a message naming the parameter - before any data is
read or written.
