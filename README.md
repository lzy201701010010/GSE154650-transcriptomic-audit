# GSE154650 transcriptomic sensitivity audit

Version: `v1.0.0`
Release date: `2026-07-20`

## Purpose

This reproducibility package supports a patient-level sensitivity reanalysis of public GSE154650 Combo-Seq data. It preserves the frozen inference boundary: disease-effect estimates are sensitive to an expression-derived histological-zone classification. It does not claim a biomarker, regulatory network, diagnostic signature, disease mechanism, causal mechanism or independently measured anatomical-zone effect.

## Data source and cohort

Raw and processed source data are controlled by NCBI GEO under `GSE154650`, with raw sequencing linked through SRA `SRP272387`. The statistical unit is the patient/library: 20 hemorrhoidal-disease cases and 18 **anal fissure surgical controls**. Each total-RNA Combo-Seq library yielded paired mRNA and miRNA measurements.

This is a secondary analysis of public, de-identified data. The repository does not redistribute patient-level raw sequencing data or the complete processed per-patient source-file set.

The deposited histological-zone variable was inferred from mRNA expression signatures, single-sample gene-set enrichment analysis and clustering. It is an **expression-derived histological-zone classification**, not an independent pathological or anatomical assessment.

## Scientific claim boundary

The package supports an inference-sensitivity conclusion only. It does not establish an anatomical cause, prove that zone explains all no-zone differences, exclude genuine low-precision disease effects, identify a biomarker, integrate mRNA with miRNA, or infer a regulatory network.

## Repository contents

- `scripts/`: path-parameterized formal analysis and figure-assembly workflows.
- `config/`: external-input requirements.
- `environment/`: frozen software and package-version records.
- `metadata/`: sanitized cohort metadata, provenance, manifests and SHA-256 indexes.
- `figure_data/`: frozen source data used for the main figures.
- `figures/`: final publication figures in SVG, PDF and TIFF formats.
- `tables/`: final Markdown tables plus machine-readable source tables.
- `supplementary/`: frozen supplementary figures and tables.
- `docs/`: final figure legends, a citation-complete manuscript copy and reference files.

## Software environment

The frozen run used R 4.4.2 on Windows 11 x64. Exact package versions are in `environment/package_versions.csv`, and complete session records are in `environment/sessionInfo.txt`. The local run was not validated across operating systems.

## Obtain the external inputs

This repository intentionally omits the GEO raw archive, the 38 processed per-patient source files and Ensembl/miRBase annotation downloads. Obtain the data from GEO `GSE154650` and linked SRA `SRP272387`. Place the 38 verified source files under `inputs/`, using the filenames in `HD_DATASET_AND_SAMPLE_MANIFEST.csv`. Supply Ensembl GRCh38 release 100 and miRBase 22.1 audit files under `source_metadata/HD_T1_annotation_cache/`. Confirm every expected SHA-256 before execution; see `config/EXTERNAL_INPUTS.md`.

## Run order

Run from the repository root with `HD_PROJECT_ROOT` set to the repository root, or make the repository root the current working directory:

1. `Rscript scripts/01_import_qc.R .`
2. `Rscript scripts/02_mrna_sensitivity_analysis.R`
3. `Rscript scripts/03_mirna_sensitivity_analysis.R`
4. `Rscript scripts/04_assemble_figures.R .`

The workflows prepare audited matrices, run the frozen mRNA and miRNA sensitivity analyses, and assemble figures from frozen outputs. No stage automatically downloads data. The public package was assembled from already frozen outputs; these scripts were not rerun during release preparation.

## Expected outputs

Formal runs create corresponding analysis-result directories. The approved publication assets are retained separately in `figures/`, `figure_data/`, `tables/` and `supplementary/`; scripts must not overwrite those package assets.

## Integrity verification

Use `python scripts/verify_sha256.py metadata/input_sha256.csv` before analysis and `python scripts/verify_sha256.py metadata/frozen_output_sha256.csv` for packaged frozen outputs. `metadata/repository_manifest.csv` and `metadata/repository_sha256.csv` inventory the complete public candidate, excluding the two self-referential index files themselves.

## Key model definition

The primary full model is `expression ~ disease group + sex + centered BMI + expression-derived histological-zone classification`. The isolated classification-addition comparison is the sex- and BMI-adjusted model without zone versus this full model. Age and batch were unavailable and were not imputed.

## Known limitations

This is a single public cohort of 38 patients. Controls are anal fissure surgical controls, not healthy or normal tissue. Zone overlap is sparse and asymmetric. The zone variable is mRNA-derived. Age and batch were unavailable. The paired miRNA layer is not an independent replication. A zero adjusted FDR count does not prove biological equivalence.

## Licence, citation and contact

Original code is released under the MIT License. GSE154650 data and third-party dependencies retain their own terms; see `LICENSE_SCOPE.md`. Citation metadata for the four confirmed authors are in `CITATION.cff`. Correspondence: Yunxiang Wu (`Wyx841106@163.com`).
