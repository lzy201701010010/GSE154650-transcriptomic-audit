#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

project_root <- normalizePath(
  Sys.getenv("HD_PROJECT_ROOT", getwd()),
  winslash = "/",
  mustWork = TRUE
)
local_lib <- file.path(project_root, ".Rlib")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(digest)
})

set.seed(20260716)

result_root <- file.path(project_root, "results", "HD_T3")
dirs <- c(
  "00_logs", "01_inputs", "02_filtering", "03_models", "04_sensitivity",
  "05_stability", "06_tables", "07_figures", "08_environment", "09_gate"
)
for (d in dirs) dir.create(file.path(result_root, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(result_root, "00_logs", "HD_T3_RUN.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
write_md <- function(path, lines) writeLines(enc2utf8(lines), path, useBytes = TRUE)
write_csv <- function(x, path) fwrite(x, path, quote = TRUE, na = "NA")
write_gz_csv <- function(x, path) fwrite(x, path, quote = TRUE, na = "NA", compress = "gzip")
sha256 <- function(path) toupper(digest(path, algo = "sha256", file = TRUE))
safe_cor <- function(x, y, method = "spearman") {
  if (sum(is.finite(x) & is.finite(y)) < 3L) return(NA_real_)
  suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
}

run_start <- Sys.time()
log_msg("HD-T3 start; exploratory-only hard stop active; HD-T4 autostart prohibited.")

paths <- list(
  manifest = file.path(project_root, "HD_DATASET_AND_SAMPLE_MANIFEST.csv"),
  raw_mirna = file.path(project_root, "results", "HD_T1", "02_matrices", "HD_T1_RAW_MIRNA.tsv.gz"),
  rpm_mirna = file.path(project_root, "results", "HD_T1", "02_matrices", "HD_T1_RPM_MIRNA.tsv.gz"),
  name_map = file.path(project_root, "results", "HD_T1", "03_feature_annotation", "HD_T1_MIRNA_NAME_MAP.csv"),
  t1_qc = file.path(project_root, "results", "HD_T1", "04_qc_tables", "HD_T1_SAMPLE_QC_METRICS.csv"),
  t1_warnings = file.path(project_root, "results", "HD_T1", "04_qc_tables", "HD_T1_SAMPLE_QC_WARNINGS.csv"),
  t1_design = file.path(project_root, "results", "HD_T1", "HD_T1_DESIGN_MATRIX_AUDIT.md"),
  t1_gate = file.path(project_root, "results", "HD_T1", "HD_T1_FINAL_GATE.md"),
  t1_output_index = file.path(project_root, "results", "HD_T1", "HD_T1_OUTPUT_INDEX.csv"),
  t1_integrity = file.path(project_root, "results", "HD_T1", "HD_T1_INPUT_INTEGRITY_REPORT.md"),
  t2_gate = file.path(project_root, "results", "HD_T2", "HD_T2_FINAL_GATE.md"),
  r2_decision = file.path(project_root, "HD_R2_FINAL_DECISION.md"),
  r2_amendment = file.path(project_root, "HD_R2_HD_T3_PROTOCOL_AMENDMENT.md"),
  prompt_review = file.path(project_root, "HD_T3_PROMPT_REVIEW.md")
)
if (!all(file.exists(unlist(paths)))) {
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  stop("Missing required HD-T3 inputs: ", paste(missing, collapse = "; "))
}

input_hashes <- data.table(
  input_name = names(paths),
  path = normalizePath(unlist(paths), winslash = "/", mustWork = TRUE),
  file_size_bytes = as.numeric(file.info(unlist(paths))$size),
  sha256 = vapply(unlist(paths), sha256, character(1))
)
write_csv(input_hashes, file.path(result_root, "01_inputs", "HD_T3_INPUT_SHA256.csv"))

manifest <- fread(paths$manifest)
raw_dt <- fread(paths$raw_mirna)
rpm_dt <- fread(paths$rpm_mirna)
name_map <- fread(paths$name_map)
qc <- fread(paths$t1_qc)
warnings <- fread(paths$t1_warnings)
t1_index <- fread(paths$t1_output_index)

required_manifest_cols <- c("patient_id", "GEO_sample_id", "disease_group", "sex", "BMI", "tissue_zone")
if (!all(required_manifest_cols %in% names(manifest))) stop("Manifest columns are incomplete.")
meta <- manifest[, ..required_manifest_cols]
meta[, disease_group := factor(
  disease_group,
  levels = c("anal fissure surgical control", "hemorrhoidal disease")
)]
meta[, sex := factor(sex, levels = c("female", "male"))]
meta[, tissue_zone := factor(
  tissue_zone,
  levels = c("anoderm region", "transition zone", "intestinal mucosa region")
)]
meta[, BMI := as.numeric(BMI)]
meta[, BMI_centered := BMI - mean(BMI)]

patient_cols <- setdiff(names(raw_dt), "feature_id")
raw_mat <- as.matrix(raw_dt[, ..patient_cols])
rpm_mat <- as.matrix(rpm_dt[, ..patient_cols])
storage.mode(raw_mat) <- "double"
storage.mode(rpm_mat) <- "double"
rownames(raw_mat) <- raw_dt$feature_id
rownames(rpm_mat) <- rpm_dt$feature_id

name_map <- name_map[match(raw_dt$feature_id, original_miRNA_id)]
if (!identical(name_map$original_miRNA_id, raw_dt$feature_id)) {
  stop("miRNA annotation order does not match raw matrix.")
}

name_map[, feature_family := fcase(
  mapping_status == "RECOGNIZED_MIRBASE_22_1_MATURE", "MATURE",
  mapping_status == "RECOGNIZED_MIRBASE_22_1_PRECURSOR", "PRECURSOR",
  default = "UNRESOLVED"
)]
name_map[, miRNA_type := fcase(
  feature_family == "MATURE", "mature_miRNA",
  feature_family == "PRECURSOR", "precursor_miRNA",
  default = "deprecated_or_unresolved_miRNA"
)]
name_map[, arm := fifelse(is.na(arm_5p_3p) | !nzchar(arm_5p_3p), NA_character_, arm_5p_3p)]
name_map[, legacy_name_status := fcase(
  mapping_status != "UNRECOGNIZED_IN_MIRBASE_22_1", "CURRENT_MIRBASE_22_1",
  (!is.na(deprecated_name) & nzchar(deprecated_name)) |
    (!is.na(deprecated_id) & nzchar(deprecated_id)), "DEPRECATED_OR_LEGACY",
  default = "UNRESOLVED_LEGACY_OR_UNVERSIONED"
)]

expected_hash <- function(path) {
  idx <- t1_index[tolower(gsub("\\\\", "/", output_path)) == tolower(gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = TRUE)))]
  if (nrow(idx) != 1L) return(NA_character_)
  toupper(idx$sha256[[1]])
}
manifest_integrity_text <- paste(readLines(paths$t1_integrity, warn = FALSE), collapse = "\n")
manifest_expected <- regmatches(
  manifest_integrity_text,
  regexpr("[A-F0-9]{64}", manifest_integrity_text)
)

integrity <- list(
  patients_38 = nrow(meta) == 38L,
  groups_20_18 = identical(as.integer(table(meta$disease_group)), c(18L, 20L)),
  raw_dimensions = nrow(raw_dt) == 1806L && ncol(raw_dt) == 39L,
  rpm_dimensions = nrow(rpm_dt) == 1806L && ncol(rpm_dt) == 39L,
  columns_match_manifest = identical(patient_cols, meta$patient_id),
  rpm_columns_match_raw = identical(names(rpm_dt), names(raw_dt)),
  raw_feature_order_matches_map = identical(raw_dt$feature_id, name_map$original_miRNA_id),
  unique_original_ids = !anyDuplicated(raw_dt$feature_id),
  unique_normalized_ids = !anyDuplicated(name_map$normalized_miRNA_id),
  no_duplicate_patients = !anyDuplicated(meta$patient_id),
  no_missing_covariates = !anyNA(meta[, .(disease_group, sex, BMI, tissue_zone)]),
  raw_no_na = !anyNA(raw_mat),
  raw_no_inf = all(is.finite(raw_mat)),
  raw_no_negative = all(raw_mat >= 0),
  rpm_no_na = !anyNA(rpm_mat),
  rpm_no_inf = all(is.finite(rpm_mat)),
  rpm_no_negative = all(rpm_mat >= 0),
  fractional_precision_present = any(abs(raw_mat - round(raw_mat)) > 1e-8),
  recognized_1776 = sum(name_map$feature_family != "UNRESOLVED") == 1776L,
  lowercase_hsa_mir_precursor_81 = sum(
    name_map$feature_family == "PRECURSOR" &
      grepl("^hsa-mir-", name_map$normalized_miRNA_id)
  ) == 81L,
  gsm4676457_mirna_not_zero_filled = {
    j <- match(meta[GEO_sample_id == "GSM4676457", patient_id], colnames(raw_mat))
    length(j) == 1L && !is.na(j) && sum(raw_mat[, j] > 0) > 0L
  },
  raw_hash_matches_t1 = identical(sha256(paths$raw_mirna), expected_hash(paths$raw_mirna)),
  rpm_hash_matches_t1 = identical(sha256(paths$rpm_mirna), expected_hash(paths$rpm_mirna)),
  map_hash_matches_t1 = identical(sha256(paths$name_map), expected_hash(paths$name_map)),
  manifest_hash_matches_t1 = length(manifest_expected) == 1L &&
    identical(sha256(paths$manifest), toupper(manifest_expected))
)
integrity_dt <- data.table(
  check = names(integrity),
  pass = vapply(integrity, isTRUE, logical(1))
)
write_csv(integrity_dt, file.path(result_root, "01_inputs", "HD_T3_INPUT_INTEGRITY_CHECKS.csv"))
if (!all(integrity_dt$pass)) {
  stop("HD_T3_FAIL_INPUT_OR_ID: ", paste(integrity_dt[pass == FALSE, check], collapse = ", "))
}
log_msg("Input integrity passed: 1806 miRNAs x 38 patients; 20 cases and 18 controls.")

design_primary <- model.matrix(~ disease_group + sex + BMI_centered + tissue_zone, data = meta)
design_unadjusted <- model.matrix(~ disease_group, data = meta)
design_no_zone <- model.matrix(~ disease_group + sex + BMI_centered, data = meta)
if (qr(design_primary)$rank != ncol(design_primary)) stop("HD_T3_FAIL_MODEL: primary design is not full rank.")
write_csv(
  data.table(patient_id = meta$patient_id, as.data.table(design_primary)),
  file.path(result_root, "01_inputs", "HD_T3_PRIMARY_DESIGN_MATRIX.csv")
)

families <- c("MATURE", "PRECURSOR", "UNRESOLVED")
filter_rows <- list()
counts_primary <- list()
counts_relaxed <- list()
counts_stringent <- list()
filter_summary_dt <- list()
for (fam in families) {
  idx <- which(name_map$feature_family == fam)
  fam_counts <- raw_mat[idx, , drop = FALSE]
  dge <- DGEList(counts = fam_counts)
  keep_primary <- filterByExpr(dge, design = design_primary)
  keep_relaxed <- filterByExpr(
    dge, design = design_primary, min.count = 5, min.total.count = 10
  )
  keep_stringent <- filterByExpr(
    dge, design = design_primary, min.count = 20, min.total.count = 30
  )
  counts_primary[[fam]] <- fam_counts[keep_primary, , drop = FALSE]
  counts_relaxed[[fam]] <- fam_counts[keep_relaxed, , drop = FALSE]
  counts_stringent[[fam]] <- fam_counts[keep_stringent, , drop = FALSE]
  audit <- copy(name_map[idx])
  audit[, `:=`(
    total_raw = rowSums(fam_counts),
    nonzero_patient_n = rowSums(fam_counts > 0),
    keep_primary_filterByExpr = keep_primary,
    keep_relaxed_min5_total10 = keep_relaxed,
    keep_stringent_min20_total30 = keep_stringent,
    primary_filter_reason = ifelse(
      keep_primary, "RETAINED_PRIMARY_FILTER",
      "LOW_EXPRESSION_FILTERBYEXPR_PRIMARY"
    )
  )]
  filter_rows[[fam]] <- audit
  filter_summary_dt[[fam]] <- data.table(
    feature_family = fam,
    before_n = length(idx),
    primary_after_n = sum(keep_primary),
    relaxed_after_n = sum(keep_relaxed),
    stringent_after_n = sum(keep_stringent)
  )
}
filter_audit <- rbindlist(filter_rows, use.names = TRUE, fill = TRUE)
filter_audit[, source_order := match(original_miRNA_id, raw_dt$feature_id)]
setorder(filter_audit, source_order)
filter_audit[, source_order := NULL]
write_csv(filter_audit, file.path(result_root, "HD_T3_MIRNA_FILTER_AUDIT.csv"))
filter_counts <- rbindlist(filter_summary_dt)
write_csv(filter_counts, file.path(result_root, "02_filtering", "HD_T3_FILTER_COUNTS_BY_FAMILY.csv"))

sample_filter_metrics <- data.table(
  patient_id = meta$patient_id,
  all_miRNA_raw_total_before = colSums(raw_mat),
  all_miRNA_nonzero_before = colSums(raw_mat > 0)
)
for (fam in families) {
  sample_filter_metrics[, paste0(tolower(fam), "_raw_total_after") :=
    colSums(counts_primary[[fam]])]
  sample_filter_metrics[, paste0(tolower(fam), "_nonzero_after") :=
    colSums(counts_primary[[fam]] > 0)]
}
write_csv(sample_filter_metrics, file.path(result_root, "02_filtering", "HD_T3_SAMPLE_FILTERING_METRICS.csv"))
write_md(
  file.path(result_root, "HD_T3_FILTERING_SUMMARY.md"),
  c(
    "# HD-T3 filtering summary",
    "",
    sprintf("- miRNAs before filtering: %d.", nrow(raw_mat)),
    sprintf("- miRNAs after primary family-separated filtering: %d.", sum(filter_counts$primary_after_n)),
    sprintf("- Mature retained: %d/%d.", filter_counts[feature_family == "MATURE", primary_after_n], filter_counts[feature_family == "MATURE", before_n]),
    sprintf("- Precursor retained: %d/%d.", filter_counts[feature_family == "PRECURSOR", primary_after_n], filter_counts[feature_family == "PRECURSOR", before_n]),
    sprintf("- Unresolved/legacy retained: %d/%d.", filter_counts[feature_family == "UNRESOLVED", primary_after_n], filter_counts[feature_family == "UNRESOLVED", before_n]),
    "- Primary filter: edgeR::filterByExpr package defaults with the frozen adjusted design.",
    "- Relaxed sensitivity: min.count=5 and min.total.count=10.",
    "- Stringent sensitivity: min.count=20 and min.total.count=30.",
    "- Filtering, TMM normalization and multiplicity correction are separate by feature family.",
    "- No feature was retained because of biological interest or observed differential results."
  )
)
log_msg("Family-separated filtering completed: ", sum(filter_counts$primary_after_n), "/", nrow(raw_mat), " retained.")

annotation_for <- function(ids) {
  ann <- name_map[match(ids, original_miRNA_id)]
  ann[, .(
    original_miRNA_id,
    normalized_miRNA_id,
    miRNA_type,
    mature_or_precursor,
    arm,
    mapping_status,
    legacy_name_status,
    feature_family
  )]
}

coef_name <- "disease_grouphemorrhoidal disease"
fit_one_family <- function(counts, model_meta, formula, model_name, family, coef_pattern = coef_name, save_plot = FALSE) {
  design <- model.matrix(formula, data = model_meta)
  rank <- qr(design)$rank
  coef_candidates <- if (length(coef_pattern) == 1L && coef_pattern %in% colnames(design)) {
    coef_pattern
  } else {
    grep(coef_pattern, colnames(design), value = TRUE)
  }
  if (nrow(counts) < 2L || rank != ncol(design) || length(coef_candidates) != 1L) {
    return(list(
      status = "NOT_ESTIMABLE", family = family, model_name = model_name,
      design = design, rank = rank, result = NULL
    ))
  }
  coef_index <- match(coef_candidates, colnames(design))
  y <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  v <- voomWithQualityWeights(y, design = design, plot = FALSE, save.plot = save_plot)
  fit <- eBayes(lmFit(v, design), robust = FALSE)
  tt <- topTable(fit, coef = coef_index, number = Inf, sort.by = "none", confint = 0.95)
  sw <- v$targets$sample.weights
  if (is.null(sw)) {
    sw <- apply(v$weights, 2, median, na.rm = TRUE)
    sw <- sw / exp(mean(log(sw)))
  }
  names(sw) <- colnames(counts)
  kish_n <- (sum(sw)^2) / sum(sw^2)
  res <- data.table(
    original_miRNA_id = rownames(tt),
    logFC = tt$logFC,
    average_expression = tt$AveExpr,
    standard_error = fifelse(is.finite(tt$t) & abs(tt$t) > 0, abs(tt$logFC / tt$t), NA_real_),
    CI_95_low = tt$CI.L,
    CI_95_high = tt$CI.R,
    moderated_t = tt$t,
    P_value = tt$P.Value,
    BH_FDR = tt$adj.P.Val,
    B_statistic = tt$B,
    model_name = model_name,
    feature_family = family,
    n_total = nrow(model_meta),
    n_case = sum(model_meta$disease_group == "hemorrhoidal disease"),
    n_control = sum(model_meta$disease_group == "anal fissure surgical control"),
    effective_sample_size_kish = kish_n,
    interpretation_level = fcase(
      family == "MATURE", "PRIMARY_FORMAL_FAMILY",
      family == "PRECURSOR", "SECONDARY_FORMAL_FAMILY",
      default = "REPORT_ONLY_NOT_HEADLINE"
    )
  )
  ann <- annotation_for(res$original_miRNA_id)
  res <- merge(ann, res, by = c("original_miRNA_id", "feature_family"), all.y = TRUE, sort = FALSE)
  setcolorder(res, c(
    "original_miRNA_id", "normalized_miRNA_id", "miRNA_type",
    "mature_or_precursor", "arm", "mapping_status", "legacy_name_status",
    "feature_family", setdiff(names(res), c(
      "original_miRNA_id", "normalized_miRNA_id", "miRNA_type",
      "mature_or_precursor", "arm", "mapping_status", "legacy_name_status",
      "feature_family"
    ))
  ))
  list(
    status = "ESTIMABLE", family = family, model_name = model_name,
    design = design, rank = rank, result = res, v = v, fit = fit,
    coef_index = coef_index, sample_weights = sw,
    norm_factors = y$samples$norm.factors
  )
}

fit_bundle <- function(count_list, model_meta, formula, model_name, coef_pattern = coef_name, save_mature_plot = FALSE) {
  fam_models <- lapply(families, function(fam) {
    fit_one_family(
      count_list[[fam]], model_meta, formula, model_name, fam,
      coef_pattern = coef_pattern,
      save_plot = save_mature_plot && fam == "MATURE"
    )
  })
  names(fam_models) <- families
  results <- rbindlist(lapply(fam_models, function(x) x$result), fill = TRUE)
  list(
    status = if (fam_models$MATURE$status == "ESTIMABLE") "ESTIMABLE" else "PRIMARY_FAMILY_NOT_ESTIMABLE",
    model_name = model_name,
    family_models = fam_models,
    result = results
  )
}

primary <- fit_bundle(
  counts_primary, meta,
  ~ disease_group + sex + BMI_centered + tissue_zone,
  "PRIMARY_ADJUSTED", save_mature_plot = TRUE
)
if (primary$status != "ESTIMABLE") stop("HD_T3_FAIL_MODEL: mature primary family is not estimable.")
unadjusted <- fit_bundle(counts_primary, meta, ~ disease_group, "UNADJUSTED")
no_zone <- fit_bundle(
  counts_primary, meta,
  ~ disease_group + sex + BMI_centered,
  "NO_ZONE_ADJUSTMENT"
)
reduced_idx <- meta$tissue_zone != "anoderm region"
reduced_counts <- lapply(counts_primary, function(x) x[, reduced_idx, drop = FALSE])
reduced <- fit_bundle(
  reduced_counts, droplevels(meta[reduced_idx]),
  ~ disease_group + sex + BMI_centered + tissue_zone,
  "REDUCED_ZONE_EXCLUDING_ANODERM"
)
transition_idx <- meta$tissue_zone == "transition zone"
transition_counts <- lapply(counts_primary, function(x) x[, transition_idx, drop = FALSE])
transition <- fit_bundle(
  transition_counts, droplevels(meta[transition_idx]),
  ~ disease_group + sex + BMI_centered,
  "TRANSITION_ZONE_LOW_PRECISION"
)
intestinal_idx <- meta$tissue_zone == "intestinal mucosa region"
intestinal_counts <- lapply(counts_primary, function(x) x[, intestinal_idx, drop = FALSE])
intestinal <- fit_bundle(
  intestinal_counts, droplevels(meta[intestinal_idx]),
  ~ disease_group + sex + BMI_centered,
  "INTESTINAL_MUCOSA_ZONE_LOW_PRECISION"
)

models <- list(
  primary = primary, unadjusted = unadjusted, no_zone = no_zone,
  reduced = reduced, transition = transition, intestinal = intestinal
)
design_audit <- rbindlist(lapply(models, function(bundle) {
  rbindlist(lapply(bundle$family_models, function(obj) data.table(
    model_name = obj$model_name,
    feature_family = obj$family,
    status = obj$status,
    n_patients = nrow(obj$design),
    design_columns = ncol(obj$design),
    design_rank = obj$rank,
    residual_degrees_of_freedom = nrow(obj$design) - obj$rank,
    coefficient_names = paste(colnames(obj$design), collapse = ";")
  )), fill = TRUE)
}), fill = TRUE)
write_csv(design_audit, file.path(result_root, "03_models", "HD_T3_DESIGN_MATRIX_AUDIT.csv"))

formal_outputs <- list(
  HD_T3_MIRNA_DE_ADJUSTED.csv.gz = primary,
  HD_T3_MIRNA_DE_UNADJUSTED.csv.gz = unadjusted,
  HD_T3_MIRNA_DE_NO_ZONE_ADJUSTMENT.csv.gz = no_zone,
  HD_T3_MIRNA_DE_REDUCED_ZONE.csv.gz = reduced,
  HD_T3_MIRNA_DE_TRANSITION_ZONE.csv.gz = transition,
  HD_T3_MIRNA_DE_INTESTINAL_MUCOSA_ZONE.csv.gz = intestinal
)
for (nm in names(formal_outputs)) {
  write_gz_csv(formal_outputs[[nm]]$result, file.path(result_root, nm))
}
log_msg("Primary and mandatory zone models completed.")

write_md(
  file.path(result_root, "HD_T3_MODEL_SPECIFICATION.md"),
  c(
    "# HD-T3 model specification",
    "",
    "- Estimand: patient-level miRNA log2 expression difference for hemorrhoidal disease minus anal fissure surgical controls.",
    "- Primary model: `~ disease_group + sex + BMI_centered + tissue_zone`.",
    "- Reference levels: anal fissure surgical control; female; anoderm region.",
    "- Input: frozen non-negative fractional raw miRNA abundance without rounding.",
    "- Filtering: family-separated edgeR::filterByExpr with the primary adjusted design.",
    "- Normalization: TMM separately for mature, precursor and unresolved families.",
    "- Mean-variance model: limma::voomWithQualityWeights.",
    "- Inference: limma::lmFit followed by limma::eBayes(robust = FALSE).",
    "- Multiplicity: Benjamini-Hochberg separately within each feature family.",
    "- Unresolved records are report-only and cannot create a headline.",
    "- Transition-zone and intestinal-mucosa analyses are low precision; anoderm is descriptive only.",
    "- Disease-by-tissue-zone interaction was not run because HD-R2 did not explicitly authorize it.",
    "- Age and batch were excluded because they are structurally unavailable.",
    "",
    "## Frozen stability rules",
    "",
    "- Near-zero classification threshold: |log2FC| < 0.05.",
    "- Stable adjusted: family FDR <0.05; non-zero and direction-consistent across unadjusted, no-zone and reduced-zone models; direction-consistent after each source-warning exclusion and both together; reproduced in at least one estimable within-zone analysis; leave-one-out direction retention >=90%; maximum leave-one-out change <=max(0.5, |primary log2FC|).",
    "- Directionally stable low precision: all directional and patient-stability checks pass but primary family FDR is >=0.05.",
    "- Alternative filter, sex-stratified and disease-by-sex results are exploratory sensitivities and cannot replace the primary gate."
  )
)

warning_scenarios <- list(
  EXCLUDE_GSM4676447 = "J01138",
  EXCLUDE_GSM4676457 = "J01148",
  EXCLUDE_BOTH_SOURCE_WARNING_FILES = c("J01138", "J01148"),
  EXCLUDE_GSM4676440_LOW_CORRELATION = "J01131",
  EXCLUDE_GSM4676456_HIGH_PCA = "J01147",
  EXCLUDE_GSM4676460_HIGH_PCA = "J01151",
  EXCLUDE_GSM4676462_HIGH_PCA_LOW_CORRELATION = "J01153",
  EXCLUDE_GSM4676467_LOW_CORRELATION = "J01158",
  EXCLUDE_GSM4676468_HIGH_PCA_LOW_CORRELATION = "J01159"
)
warning_models <- list()
warning_long <- list()
warning_scenario_summary <- list()
for (nm in names(warning_scenarios)) {
  keep <- !meta$patient_id %chin% warning_scenarios[[nm]]
  sub_counts <- lapply(counts_primary, function(x) x[, keep, drop = FALSE])
  obj <- fit_bundle(
    sub_counts, droplevels(meta[keep]),
    ~ disease_group + sex + BMI_centered + tissue_zone,
    nm
  )
  warning_models[[nm]] <- obj
  tmp <- obj$result[, .(
    original_miRNA_id, feature_family,
    scenario = nm, scenario_logFC = logFC, scenario_FDR = BH_FDR
  )]
  primary_tmp <- primary$result[, .(
    original_miRNA_id, feature_family,
    primary_logFC = logFC, primary_FDR = BH_FDR
  )]
  tmp <- merge(primary_tmp, tmp, by = c("original_miRNA_id", "feature_family"), all.x = TRUE, sort = FALSE)
  warning_long[[nm]] <- tmp
  warning_scenario_summary[[nm]] <- data.table(
    scenario = nm,
    excluded_patient_ids = paste(warning_scenarios[[nm]], collapse = ";"),
    n_patients = sum(keep),
    mature_FDR_significant_n = sum(obj$result$feature_family == "MATURE" & obj$result$BH_FDR < 0.05, na.rm = TRUE),
    precursor_FDR_significant_n = sum(obj$result$feature_family == "PRECURSOR" & obj$result$BH_FDR < 0.05, na.rm = TRUE),
    logFC_spearman_vs_primary = safe_cor(obj$result$logFC, primary$result$logFC),
    direction_concordance_vs_primary = mean(
      sign(obj$result$logFC) == sign(primary$result$logFC), na.rm = TRUE
    ),
    median_absolute_logFC_change = median(abs(obj$result$logFC - primary$result$logFC), na.rm = TRUE)
  )
}
warning_long_dt <- rbindlist(warning_long, fill = TRUE)
warning_long_dt <- merge(
  warning_long_dt,
  annotation_for(unique(warning_long_dt$original_miRNA_id)),
  by = c("original_miRNA_id", "feature_family"), all.x = TRUE, sort = FALSE
)
write_csv(warning_long_dt, file.path(result_root, "HD_T3_WARNING_SAMPLE_SENSITIVITY.csv"))
write_csv(
  rbindlist(warning_scenario_summary),
  file.path(result_root, "04_sensitivity", "HD_T3_WARNING_SAMPLE_SCENARIO_SUMMARY.csv")
)

fit_filter_bundle <- function(count_list, nm) {
  obj <- fit_bundle(
    count_list, meta,
    ~ disease_group + sex + BMI_centered + tissue_zone,
    nm
  )
  write_gz_csv(obj$result, file.path(result_root, "04_sensitivity", paste0("HD_T3_", nm, ".csv.gz")))
  obj
}
filter_relaxed <- fit_filter_bundle(counts_relaxed, "FILTER_RELAXED_MIN5_TOTAL10")
filter_stringent <- fit_filter_bundle(counts_stringent, "FILTER_STRINGENT_MIN20_TOTAL30")

sex_models <- list()
for (sex_level in levels(meta$sex)) {
  idx <- meta$sex == sex_level
  sub_counts <- lapply(counts_primary, function(x) x[, idx, drop = FALSE])
  sex_models[[sex_level]] <- fit_bundle(
    sub_counts, droplevels(meta[idx]),
    ~ disease_group + BMI_centered + tissue_zone,
    paste0("SEX_STRATIFIED_", toupper(sex_level))
  )
  write_gz_csv(
    sex_models[[sex_level]]$result,
    file.path(result_root, "04_sensitivity", paste0("HD_T3_MIRNA_DE_SEX_", toupper(sex_level), ".csv.gz"))
  )
}

sex_interaction <- fit_bundle(
  counts_primary, meta,
  ~ disease_group * sex + BMI_centered + tissue_zone,
  "EXPLORATORY_DISEASE_BY_SEX_INTERACTION",
  coef_pattern = "disease_group.*:sex"
)
write_gz_csv(
  sex_interaction$result,
  file.path(result_root, "04_sensitivity", "HD_T3_DISEASE_BY_SEX_INTERACTION.csv.gz")
)

primary_ids <- primary$result$original_miRNA_id
primary_effect <- setNames(primary$result$logFC, primary_ids)
primary_se <- setNames(primary$result$standard_error, primary_ids)
loo_logfc <- matrix(
  NA_real_, nrow = length(primary_ids), ncol = nrow(meta),
  dimnames = list(primary_ids, meta$patient_id)
)
loo_fdr <- matrix(
  NA_real_, nrow = length(primary_ids), ncol = nrow(meta),
  dimnames = list(primary_ids, meta$patient_id)
)
sample_influence <- vector("list", nrow(meta))
log_msg("Starting leave-one-patient-out refits: ", nrow(meta), " iterations.")
for (i in seq_len(nrow(meta))) {
  keep <- seq_len(nrow(meta)) != i
  sub_counts <- lapply(counts_primary, function(x) x[, keep, drop = FALSE])
  obj <- fit_bundle(
    sub_counts, droplevels(meta[keep]),
    ~ disease_group + sex + BMI_centered + tissue_zone,
    paste0("LOO_", meta$patient_id[i])
  )
  idx <- match(obj$result$original_miRNA_id, primary_ids)
  loo_logfc[idx, i] <- obj$result$logFC
  loo_fdr[idx, i] <- obj$result$BH_FDR
  delta <- obj$result$logFC - primary_effect[obj$result$original_miRNA_id]
  std_delta <- abs(delta) / pmax(primary_se[obj$result$original_miRNA_id], 0.05)
  sample_influence[[i]] <- data.table(
    patient_id = meta$patient_id[i],
    GEO_sample_id = meta$GEO_sample_id[i],
    model_status = obj$status,
    median_absolute_logFC_change = median(abs(delta), na.rm = TRUE),
    q95_absolute_logFC_change = as.numeric(quantile(abs(delta), 0.95, na.rm = TRUE)),
    maximum_absolute_logFC_change = max(abs(delta), na.rm = TRUE),
    cook_like_max_standardized_change = max(std_delta, na.rm = TRUE),
    direction_flip_fraction = mean(sign(obj$result$logFC) != sign(primary_effect[obj$result$original_miRNA_id]), na.rm = TRUE),
    rank_spearman_vs_primary = safe_cor(
      rank(abs(obj$result$moderated_t), ties.method = "average"),
      rank(abs(primary$result$moderated_t[idx]), ties.method = "average")
    )
  )
  if (i %% 5L == 0L || i == nrow(meta)) log_msg("LOO progress ", i, "/", nrow(meta), ".")
}

near_zero <- function(x) !is.na(x) & abs(x) < 0.05
direction_agree <- function(a, b) {
  !is.na(a) & !is.na(b) & (
    (!near_zero(a) & !near_zero(b) & sign(a) == sign(b)) |
      (near_zero(a) & near_zero(b))
  )
}
strict_direction_agree <- function(a, b) {
  !is.na(a) & !is.na(b) & !near_zero(a) & !near_zero(b) & sign(a) == sign(b)
}
loo_same <- sweep(loo_logfc, 1, primary_effect[rownames(loo_logfc)], function(x, p) {
  direction_agree(x, p)
})
loo_summary <- data.table(
  original_miRNA_id = rownames(loo_logfc),
  primary_logFC = primary_effect[rownames(loo_logfc)],
  direction_retention_fraction = rowMeans(loo_same, na.rm = TRUE),
  loo_min_logFC = apply(loo_logfc, 1, min, na.rm = TRUE),
  loo_max_logFC = apply(loo_logfc, 1, max, na.rm = TRUE),
  loo_max_absolute_change = apply(
    abs(loo_logfc - primary_effect[rownames(loo_logfc)]), 1, max, na.rm = TRUE
  ),
  loo_median_absolute_change = apply(
    abs(loo_logfc - primary_effect[rownames(loo_logfc)]), 1, median, na.rm = TRUE
  ),
  loo_FDR_lt_0_05_fraction = rowMeans(loo_fdr < 0.05, na.rm = TRUE)
)
loo_summary <- merge(
  annotation_for(loo_summary$original_miRNA_id),
  loo_summary, by = "original_miRNA_id", all.y = TRUE, sort = FALSE
)
write_csv(loo_summary, file.path(result_root, "HD_T3_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv"))
write_gz_csv(
  data.table(original_miRNA_id = rownames(loo_logfc), as.data.table(loo_logfc)),
  file.path(result_root, "04_sensitivity", "HD_T3_LEAVE_ONE_PATIENT_OUT_LOGFC_MATRIX.csv.gz")
)

quality_weights <- data.table(patient_id = meta$patient_id)
for (fam in families) {
  obj <- primary$family_models[[fam]]
  qcol <- paste0(tolower(fam), "_quality_weight")
  ncol_name <- paste0(tolower(fam), "_normalization_factor")
  quality_weights[, (qcol) := NA_real_]
  quality_weights[, (ncol_name) := NA_real_]
  if (obj$status == "ESTIMABLE") {
    quality_weights[, (qcol) := as.numeric(obj$sample_weights[patient_id])]
    quality_weights[, (ncol_name) := as.numeric(obj$norm_factors)]
  }
}
quality_weights[, quality_weight := mature_quality_weight]
quality_weights <- merge(quality_weights, meta, by = "patient_id", all.x = TRUE, sort = FALSE)
quality_weights <- merge(
  quality_weights,
  warnings[, .(patient_id, warning_class, warning_reasons, primary_analysis_status)],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
quality_weights <- merge(
  quality_weights,
  qc[analyte == "miRNA", .(
    patient_id, raw_total, nonzero_feature_n, detection_fraction,
    median_spearman_to_others, minimum_spearman_to_others,
    PCA_distance, PC1, PC2
  )],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
quality_weights[, T1_miRNA_warning := grepl("miRNA_", warning_reasons, fixed = TRUE)]
write_csv(quality_weights, file.path(result_root, "HD_T3_SAMPLE_QUALITY_WEIGHTS.csv"))

sample_influence_dt <- rbindlist(sample_influence, fill = TRUE)
sample_influence_dt <- merge(
  sample_influence_dt, quality_weights,
  by = c("patient_id", "GEO_sample_id"), all.x = TRUE, sort = FALSE
)
sample_influence_dt[, leverage_primary_design := hat(design_primary)]
sample_influence_dt[, high_disease_effect_influence := (
  rank_spearman_vs_primary < 0.90 |
    direction_flip_fraction > 0.10 |
    cook_like_max_standardized_change > 4
)]
write_csv(sample_influence_dt, file.path(result_root, "HD_T3_SAMPLE_INFLUENCE_METRICS.csv"))

get_model_cols <- function(bundle, prefix) {
  bundle$result[, .(
    original_miRNA_id,
    feature_family,
    value_logFC = logFC,
    value_FDR = BH_FDR,
    value_t = moderated_t
  )][, c(
    "value_logFC", "value_FDR", "value_t"
  ) := NULL]
}
model_extract <- function(bundle, prefix) {
  x <- bundle$result[, .(
    original_miRNA_id, feature_family,
    logFC, BH_FDR, moderated_t
  )]
  setnames(
    x, c("logFC", "BH_FDR", "moderated_t"),
    paste0(prefix, c("_logFC", "_FDR", "_t"))
  )
  x
}
comparison <- Reduce(
  function(x, y) merge(x, y, by = c("original_miRNA_id", "feature_family"), all = TRUE, sort = FALSE),
  list(
    model_extract(primary, "primary"),
    model_extract(unadjusted, "unadjusted"),
    model_extract(no_zone, "no_zone"),
    model_extract(reduced, "reduced_zone"),
    model_extract(transition, "transition_zone"),
    model_extract(intestinal, "intestinal_zone"),
    loo_summary[, .(
      original_miRNA_id, feature_family,
      direction_retention_fraction, loo_max_absolute_change,
      loo_median_absolute_change, loo_FDR_lt_0_05_fraction
    )]
  )
)

warning_wide_logfc <- dcast(
  warning_long_dt,
  original_miRNA_id + feature_family ~ scenario,
  value.var = "scenario_logFC"
)
comparison <- merge(
  comparison, warning_wide_logfc,
  by = c("original_miRNA_id", "feature_family"), all.x = TRUE, sort = FALSE
)
comparison <- merge(
  annotation_for(comparison$original_miRNA_id),
  comparison, by = c("original_miRNA_id", "feature_family"), all.y = TRUE, sort = FALSE
)

source_warning_cols <- c(
  "EXCLUDE_GSM4676447",
  "EXCLUDE_GSM4676457",
  "EXCLUDE_BOTH_SOURCE_WARNING_FILES"
)
all_warning_cols <- names(warning_scenarios)
comparison[, source_warning_direction_stable := Reduce(
  `&`, lapply(.SD, function(x) strict_direction_agree(primary_logFC, x))
), .SDcols = source_warning_cols]
comparison[, all_warning_direction_stable := Reduce(
  `&`, lapply(.SD, function(x) strict_direction_agree(primary_logFC, x))
), .SDcols = all_warning_cols]
comparison[, within_zone_direction_reproduced := (
  strict_direction_agree(primary_logFC, transition_zone_logFC) |
    strict_direction_agree(primary_logFC, intestinal_zone_logFC)
)]
comparison[, core_model_direction_stable := (
  strict_direction_agree(primary_logFC, unadjusted_logFC) &
    strict_direction_agree(primary_logFC, no_zone_logFC) &
    strict_direction_agree(primary_logFC, reduced_zone_logFC)
)]
comparison[, sample_sensitive := (
  direction_retention_fraction < 0.90 |
    loo_max_absolute_change > pmax(0.5, abs(primary_logFC)) |
    !all_warning_direction_stable
)]
comparison[, zone_dependent := (
  !direction_agree(primary_logFC, no_zone_logFC) |
    !direction_agree(primary_logFC, reduced_zone_logFC) |
    (
      abs(no_zone_logFC) >= 0.10 &
        abs(primary_logFC) < 0.50 * abs(no_zone_logFC)
    )
)]
comparison[, model_dependent := (
  !direction_agree(primary_logFC, unadjusted_logFC) |
    !direction_agree(primary_logFC, no_zone_logFC)
)]
comparison[, stable_without_fdr := (
  core_model_direction_stable &
    source_warning_direction_stable &
    within_zone_direction_reproduced &
    direction_retention_fraction >= 0.90 &
    loo_max_absolute_change <= pmax(0.5, abs(primary_logFC))
)]
comparison[, adjusted_stable := primary_FDR < 0.05 & stable_without_fdr]
comparison[, low_precision_stable := primary_FDR >= 0.05 & stable_without_fdr]
comparison[, only_unadjusted := (
  unadjusted_FDR < 0.05 & primary_FDR >= 0.05
)]

comparison[, stability_class := fcase(
  adjusted_stable, "ADJUSTED_STABLE",
  low_precision_stable, "DIRECTIONALLY_STABLE_LOW_PRECISION",
  sample_sensitive, "SAMPLE_SENSITIVE",
  zone_dependent, "ZONE_DEPENDENT",
  only_unadjusted, "ONLY_UNADJUSTED",
  model_dependent, "MODEL_DEPENDENT",
  default = "UNSTABLE"
)]

filtered_out <- filter_audit[keep_primary_filterByExpr == FALSE]
filtered_class <- annotation_for(filtered_out$original_miRNA_id)
filtered_class[, stability_class := "INSUFFICIENT_EXPRESSION"]
stability_out <- rbindlist(list(comparison, filtered_class), fill = TRUE)
stability_out[, source_order := match(original_miRNA_id, raw_dt$feature_id)]
setorder(stability_out, source_order)
stability_out[, source_order := NULL]
write_csv(stability_out, file.path(result_root, "HD_T3_MIRNA_STABILITY_CLASSIFICATION.csv"))

comparison[, primary_rank := frank(primary_FDR, ties.method = "average")]
comparison[, unadjusted_rank := frank(unadjusted_FDR, ties.method = "average")]
comparison[, no_zone_rank := frank(no_zone_FDR, ties.method = "average")]
comparison[, adjusted_minus_unadjusted_logFC := primary_logFC - unadjusted_logFC]
comparison[, adjusted_minus_no_zone_logFC := primary_logFC - no_zone_logFC]
comparison[, direction_reversed_after_adjustment := !direction_agree(primary_logFC, no_zone_logFC)]
comparison[, zone_effect_retention_ratio := fifelse(
  abs(no_zone_logFC) > 0.05, abs(primary_logFC) / abs(no_zone_logFC), NA_real_
)]
comparison[, zone_dependence_class := fcase(
  is.na(primary_logFC) | is.na(no_zone_logFC), "UNINTERPRETABLE",
  !direction_agree(primary_logFC, no_zone_logFC) |
    !direction_agree(primary_logFC, reduced_zone_logFC), "DIRECTION_REVERSED",
  abs(no_zone_logFC) >= 0.10 & abs(primary_logFC) < 0.50 * abs(no_zone_logFC), "LARGELY_ZONE_EXPLAINED",
  abs(no_zone_logFC) >= 0.10 & abs(primary_logFC) < 0.80 * abs(no_zone_logFC), "PARTIALLY_ATTENUATED_AFTER_ZONE_ADJUSTMENT",
  strict_direction_agree(primary_logFC, no_zone_logFC) &
    strict_direction_agree(primary_logFC, reduced_zone_logFC) &
    abs(primary_logFC) >= 0.80 * abs(no_zone_logFC), "RELATIVELY_ZONE_INDEPENDENT",
  direction_agree(primary_logFC, no_zone_logFC), "DIRECTION_PRESERVED_LOW_PRECISION",
  default = "UNINTERPRETABLE"
)]
write_csv(comparison, file.path(result_root, "HD_T3_ZONE_DEPENDENCE_AUDIT.csv"))
log_msg("Stability and tissue-zone dependence classifications completed.")

theme_hd <- theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    legend.title = element_text(size = 6.7),
    legend.text = element_text(size = 6.2),
    strip.text = element_text(size = 6.7, face = "bold"),
    plot.title = element_text(size = 7.5, face = "bold"),
    panel.grid = element_blank()
  )
palette_group <- c(
  "anal fissure surgical control" = "#3182BD",
  "hemorrhoidal disease" = "#D24B40"
)
save_plot_pair <- function(plot, id, width_mm = 165, height_mm = 115) {
  pdf_path <- file.path(result_root, "07_figures", paste0(id, ".pdf"))
  png_path <- file.path(result_root, "07_figures", paste0(id, ".png"))
  ggsave(pdf_path, plot, width = width_mm / 25.4, height = height_mm / 25.4, device = cairo_pdf)
  ggsave(png_path, plot, width = width_mm / 25.4, height = height_mm / 25.4, dpi = 600, bg = "white")
  data.table(
    figure_id = id, pdf_path = pdf_path, png_path = png_path,
    width_mm = width_mm, height_mm = height_mm, dpi_png = 600
  )
}
figure_index <- list()

p1 <- ggplot(filter_counts, aes(feature_family, before_n, fill = "Before")) +
  geom_col(width = 0.62, fill = "#D8D8D8") +
  geom_col(aes(y = primary_after_n, fill = "After"), width = 0.42) +
  scale_fill_manual(values = c(Before = "#D8D8D8", After = "#3182BD")) +
  labs(title = "Family-separated miRNA filtering", x = NULL, y = "miRNA features", fill = NULL) +
  theme_hd
figure_index[[1]] <- save_plot_pair(p1, "HD_T3_FIG01_FILTER_COUNTS", 120, 90)

v_mature <- primary$family_models$MATURE$v
if (!is.null(v_mature$voom.xy) && !is.null(v_mature$voom.line)) {
  voom_dt <- data.table(
    mean_log_count = v_mature$voom.xy$x,
    sqrt_residual_sd = v_mature$voom.xy$y
  )
  line_dt <- data.table(
    mean_log_count = v_mature$voom.line$x,
    sqrt_residual_sd = v_mature$voom.line$y
  )
  p2 <- ggplot(voom_dt, aes(mean_log_count, sqrt_residual_sd)) +
    geom_point(size = 0.55, alpha = 0.30, colour = "#767676") +
    geom_line(data = line_dt, colour = "#D24B40", linewidth = 0.65) +
    labs(title = "Mature-miRNA voom mean-variance trend", x = "Mean log2 count", y = "Square-root residual SD") +
    theme_hd
  figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p2, "HD_T3_FIG02_VOOM_MEAN_VARIANCE")
}

qw_plot <- copy(quality_weights)[order(quality_weight)]
qw_plot[, patient_id := factor(patient_id, levels = patient_id)]
p3 <- ggplot(qw_plot, aes(patient_id, quality_weight, fill = disease_group)) +
  geom_col(width = 0.78) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
  scale_fill_manual(values = palette_group) +
  labs(title = "Mature-miRNA sample quality weights", x = "Patient", y = "Quality weight", fill = NULL) +
  theme_hd + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p3, "HD_T3_FIG03_SAMPLE_QUALITY_WEIGHTS", 183, 100)

p4 <- ggplot(quality_weights, aes(PC1, PC2, colour = disease_group, shape = tissue_zone)) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_colour_manual(values = palette_group) +
  labs(title = "T1 miRNA PCA with frozen patient labels", colour = NULL, shape = "Tissue zone") +
  theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p4, "HD_T3_FIG04_MIRNA_PCA", 150, 105)

p5 <- ggplot(comparison, aes(unadjusted_logFC, primary_logFC, colour = feature_family)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#767676", linewidth = 0.4) +
  geom_point(alpha = 0.50, size = 0.75) +
  scale_colour_manual(values = c(MATURE = "#3182BD", PRECURSOR = "#E28E2C", UNRESOLVED = "#767676")) +
  labs(
    title = "Adjusted versus unadjusted miRNA effects",
    subtitle = sprintf("Spearman rho = %.3f", safe_cor(comparison$unadjusted_logFC, comparison$primary_logFC)),
    x = "Unadjusted log2 fold change", y = "Adjusted log2 fold change", colour = "Family"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p5, "HD_T3_FIG05_ADJUSTED_VS_UNADJUSTED")

p6 <- ggplot(comparison, aes(no_zone_logFC, primary_logFC, colour = zone_dependence_class)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#767676", linewidth = 0.4) +
  geom_point(alpha = 0.55, size = 0.75) +
  labs(
    title = "Effect change after tissue-zone adjustment",
    subtitle = sprintf("Spearman rho = %.3f", safe_cor(comparison$no_zone_logFC, comparison$primary_logFC)),
    x = "No-zone-adjustment log2 fold change", y = "Adjusted log2 fold change", colour = "Zone class"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p6, "HD_T3_FIG06_ZONE_ADJUSTMENT_EFFECT", 183, 115)

p7 <- ggplot(comparison, aes(reduced_zone_logFC, primary_logFC, colour = stability_class)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#767676", linewidth = 0.4) +
  geom_point(alpha = 0.55, size = 0.75) +
  labs(
    title = "Primary versus reduced-zone miRNA effects",
    x = "Reduced-zone log2 fold change", y = "Adjusted log2 fold change", colour = "Stability"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p7, "HD_T3_FIG07_PRIMARY_VS_REDUCED_ZONE", 183, 115)

zone_counts <- comparison[, .N, by = zone_dependence_class][order(-N)]
p8 <- ggplot(zone_counts, aes(reorder(zone_dependence_class, N), N, fill = zone_dependence_class)) +
  geom_col(show.legend = FALSE) + coord_flip() +
  labs(title = "Tissue-zone dependence classification", x = NULL, y = "Tested miRNAs") +
  theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p8, "HD_T3_FIG08_ZONE_DEPENDENCE_CLASSIFICATION", 150, 100)

stability_counts <- stability_out[, .N, by = stability_class][order(-N)]
p9 <- ggplot(stability_counts, aes(reorder(stability_class, N), N, fill = stability_class)) +
  geom_col(show.legend = FALSE) + coord_flip() +
  labs(title = "miRNA stability classification", x = NULL, y = "miRNA features") +
  theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p9, "HD_T3_FIG09_STABILITY_CLASSIFICATION", 150, 105)

candidate <- comparison[
  feature_family == "MATURE"
][order(-adjusted_stable, -low_precision_stable, primary_FDR, -abs(primary_logFC))]
top_ids <- head(candidate$original_miRNA_id, 6L)
if (length(top_ids) > 0L) {
  mature_E <- primary$family_models$MATURE$v$E
  top_ids <- top_ids[top_ids %in% rownames(mature_E)]
  expr_long <- data.table(
    original_miRNA_id = rep(top_ids, times = ncol(mature_E)),
    patient_id = rep(colnames(mature_E), each = length(top_ids)),
    voom_logCPM = as.vector(mature_E[top_ids, , drop = FALSE])
  )
  expr_long <- merge(expr_long, meta, by = "patient_id", all.x = TRUE, sort = FALSE)
  p10 <- ggplot(expr_long, aes(disease_group, voom_logCPM, colour = disease_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.35) +
    geom_jitter(width = 0.13, size = 0.85, alpha = 0.75) +
    facet_wrap(~ original_miRNA_id, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = palette_group) +
    scale_x_discrete(labels = c(
      "anal fissure surgical control" = "anal fissure\ncontrol",
      "hemorrhoidal disease" = "hemorrhoidal\ndisease"
    )) +
    labs(title = "Leading mature-miRNA patient-level expression", x = NULL, y = "Voom log2 CPM") +
    theme_hd + theme(legend.position = "none")
  figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p10, "HD_T3_FIG10_TOP_MATURE_MIRNA_EXPRESSION", 183, 125)

  loo_plot <- loo_summary[original_miRNA_id %chin% top_ids]
  loo_plot[, original_miRNA_id := factor(original_miRNA_id, levels = rev(top_ids))]
  p11 <- ggplot(loo_plot, aes(primary_logFC, original_miRNA_id)) +
    geom_errorbarh(aes(xmin = loo_min_logFC, xmax = loo_max_logFC), height = 0.20, colour = "#767676") +
    geom_point(colour = "#3182BD", size = 1.8) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
    labs(title = "Leave-one-patient-out effect ranges", x = "Adjusted log2 fold change", y = NULL) +
    theme_hd
  figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p11, "HD_T3_FIG11_LEAVE_ONE_OUT_RANGES", 150, 100)
}

figure_index_dt <- rbindlist(figure_index, fill = TRUE)
figure_index_dt[, `:=`(
  statistical_unit = "patient/library",
  status = "HD-T3 candidate analysis figure; not a frozen manuscript figure",
  visual_QA = "PENDING_MANUAL_INSPECTION"
)]
write_csv(figure_index_dt, file.path(result_root, "07_figures", "HD_T3_FIGURE_INDEX.csv"))
write_md(
  file.path(result_root, "HD_T3_FIGURE_INDEX.md"),
  c(
    "# HD-T3 figure index",
    "",
    "- Core conclusion: determine whether adjusted miRNA effects persist across tissue-zone and patient-influence constraints.",
    "- Evidence chain: filtering and mean-variance validity; sample quality/PCA; adjustment-induced effect change; reduced-zone consistency; stability classes; patient-level and leave-one-out views.",
    "- Archetype: quantitative robustness series.",
    "- Backend: R only.",
    "- Export: PDF plus 600-dpi PNG.",
    "- Review risk: sparse within-zone support, family-specific multiplicity, source-warning and high-influence patients, precursor/unresolved identity limits.",
    "- These are candidate HD-T3 analysis figures, not final manuscript figures.",
    "",
    paste0("- ", figure_index_dt$figure_id, ": ", basename(figure_index_dt$pdf_path), " + ", basename(figure_index_dt$png_path))
  )
)
write_md(
  file.path(result_root, "07_figures", "HD_T3_FIGURE_QA_NOTES.md"),
  c(
    "# HD-T3 figure visual QA",
    "",
    "- Export backend: R only.",
    "- PDF/PNG pair generated for every indexed figure.",
    "- PNG resolution: 600 dpi.",
    "- White background and restrained non-rainbow palette used.",
    "- Quantitative panels trace to HD-T3 CSV or model outputs.",
    "- Manual visual inspection: PENDING_MANUAL_INSPECTION."
  )
)
log_msg("Candidate figures exported as PDF plus 600-dpi PNG.")

primary_sig_mature <- sum(primary$result$feature_family == "MATURE" & primary$result$BH_FDR < 0.05, na.rm = TRUE)
primary_sig_precursor <- sum(primary$result$feature_family == "PRECURSOR" & primary$result$BH_FDR < 0.05, na.rm = TRUE)
unadjusted_sig_mature <- sum(unadjusted$result$feature_family == "MATURE" & unadjusted$result$BH_FDR < 0.05, na.rm = TRUE)
reduced_sig_mature <- sum(reduced$result$feature_family == "MATURE" & reduced$result$BH_FDR < 0.05, na.rm = TRUE)
stable_mature <- sum(comparison$feature_family == "MATURE" & comparison$adjusted_stable, na.rm = TRUE)
stable_precursor <- sum(comparison$feature_family == "PRECURSOR" & comparison$adjusted_stable, na.rm = TRUE)
stable_unresolved <- sum(comparison$feature_family == "UNRESOLVED" & comparison$adjusted_stable, na.rm = TRUE)
low_precision_mature <- sum(comparison$feature_family == "MATURE" & comparison$low_precision_stable, na.rm = TRUE)
zone_dependent_n <- sum(comparison$stability_class == "ZONE_DEPENDENT", na.rm = TRUE)
sample_sensitive_n <- sum(comparison$stability_class == "SAMPLE_SENSITIVE", na.rm = TRUE)
only_unadjusted_n <- sum(comparison$stability_class == "ONLY_UNADJUSTED", na.rm = TRUE)

final_status <- if (stable_mature > 0L) {
  "HD_T3_PASS_ADJUSTED_STABLE_MIRNA_SIGNAL_REQUIRES_WORK_REVIEW"
} else if (low_precision_mature > 0L) {
  "HD_T3_PASS_LIMITED_LOW_PRECISION_SIGNAL_REQUIRES_WORK_REVIEW"
} else if (only_unadjusted_n > 0L || (
  primary_sig_mature == 0L && unadjusted_sig_mature > 0L
)) {
  "HD_T3_STOP_ONLY_UNADJUSTED_OR_ZONE_DEPENDENT_SIGNAL"
} else if (sample_sensitive_n > 0L && any(
  comparison$feature_family == "MATURE" &
    comparison$primary_FDR < 0.05 &
    comparison$sample_sensitive,
  na.rm = TRUE
)) {
  "HD_T3_STOP_SAMPLE_SENSITIVE_SIGNAL"
} else {
  "HD_T3_STOP_NO_STABLE_MIRNA_SIGNAL"
}

hd_t4_eligibility <- if (stable_mature > 0L) {
  "ELIGIBLE_FOR_WORK_REVIEW_ONLY_NOT_AUTHORIZED"
} else {
  "NOT_ELIGIBLE_HARD_STOP"
}
independent_mirna_information <- if (stable_mature > 0L) {
  "FDR_SUPPORTED_STABLE_MATURE_MIRNA_ASSOCIATION_PRESENT"
} else if (low_precision_mature > 0L) {
  "DIRECTIONAL_LOW_PRECISION_INFORMATION_ONLY"
} else {
  "NO_STABLE_INFORMATION_INDEPENDENT_OF_MRNA_ZONE_DEPENDENCE"
}

results_ledger <- data.table(
  item = c(
    "patients_primary", "miRNA_before_filter", "miRNA_after_filter",
    "mature_before_filter", "mature_after_filter",
    "precursor_before_filter", "precursor_after_filter",
    "unresolved_before_filter", "unresolved_after_filter",
    "adjusted_mature_FDR_significant_n", "adjusted_precursor_FDR_significant_n",
    "unadjusted_mature_FDR_significant_n", "reduced_zone_mature_FDR_significant_n",
    "adjusted_stable_mature_n", "adjusted_stable_precursor_n",
    "adjusted_stable_unresolved_n", "low_precision_stable_mature_n",
    "zone_dependent_n", "sample_sensitive_n", "only_unadjusted_n",
    "mature_quality_weight_min", "mature_quality_weight_max",
    "independent_miRNA_information", "HD_T4_eligibility", "final_status"
  ),
  value = c(
    nrow(meta), nrow(raw_mat), sum(filter_counts$primary_after_n),
    filter_counts[feature_family == "MATURE", before_n],
    filter_counts[feature_family == "MATURE", primary_after_n],
    filter_counts[feature_family == "PRECURSOR", before_n],
    filter_counts[feature_family == "PRECURSOR", primary_after_n],
    filter_counts[feature_family == "UNRESOLVED", before_n],
    filter_counts[feature_family == "UNRESOLVED", primary_after_n],
    primary_sig_mature, primary_sig_precursor,
    unadjusted_sig_mature, reduced_sig_mature,
    stable_mature, stable_precursor, stable_unresolved,
    low_precision_mature, zone_dependent_n, sample_sensitive_n,
    only_unadjusted_n, min(quality_weights$quality_weight),
    max(quality_weights$quality_weight), independent_mirna_information,
    hd_t4_eligibility, final_status
  )
)
write_csv(results_ledger, file.path(result_root, "06_tables", "HD_T3_RESULTS_LEDGER.csv"))
write_md(
  file.path(result_root, "HD_T3_RESULTS_LEDGER.md"),
  c("# HD-T3 results ledger", "", paste0("- ", results_ledger$item, ": ", results_ledger$value))
)

write_md(
  file.path(result_root, "HD_T3_CLAIM_BOUNDARY.md"),
  c(
    "# HD-T3 claim boundary",
    "",
    "## Allowed",
    "",
    "- Patient-level adjusted miRNA associations relative to anal fissure surgical controls.",
    "- Tissue-zone sensitivity, direction stability, uncertainty and patient influence.",
    "- Mature, precursor and unresolved identity with explicit family limitations.",
    "- A transparently negative or low-precision orthogonal miRNA result.",
    "",
    "## Prohibited",
    "",
    "- miRNA causation, direct mRNA regulation, diagnostic biomarker or treatment-target claims.",
    "- Repackaging the 6,111 unadjusted mRNA genes or rescuing the HD-T2 negative feature-level result.",
    "- Target-database enrichment, miRNA-mRNA correlation/network, hub selection or single-cell mapping.",
    "- Claiming tissue-zone confounding has been eliminated.",
    "",
    paste0("HD-T3 final state: `", final_status, "`."),
    paste0("Independent miRNA information: `", independent_mirna_information, "`."),
    paste0("HD-T4 eligibility: `", hd_t4_eligibility, "`.")
  )
)

write_md(
  file.path(result_root, "HD_T3_HD_T4_ELIGIBILITY.md"),
  c(
    "# HD-T3 to HD-T4 eligibility",
    "",
    paste0("**Decision:** `", hd_t4_eligibility, "`"),
    "",
    sprintf("- Stable adjusted mature miRNAs: %d.", stable_mature),
    sprintf("- Stable adjusted precursor miRNAs: %d.", stable_precursor),
    sprintf("- Directionally stable low-precision mature miRNAs: %d.", low_precision_mature),
    paste0("- Independent miRNA information: `", independent_mirna_information, "`."),
    "- Even an eligible state requires separate Work review and explicit user authorization.",
    "- No HD-T4 directory, target database, miRNA-mRNA association or network was created."
  )
)

write_md(
  file.path(result_root, "HD_T3_INPUT_INTEGRITY_REPORT.md"),
  c(
    "# HD-T3 input integrity report",
    "",
    "- Formal raw and RPM miRNA matrices: 1,806 features x 38 patients.",
    "- Patients: 20 hemorrhoidal disease and 18 anal fissure surgical controls.",
    "- Column order exactly matches the formal manifest.",
    "- No NA, Inf, negative values, duplicate patients or duplicate original/normalized miRNA identifiers.",
    "- Fractional raw precision was preserved; no integer coercion was used.",
    "- GSM4676457 common miRNA values were not imputed or zero-filled.",
    "- Recognized miRBase 22.1 records: 1,776.",
    "- Lowercase hsa-mir precursor records: 81; all precursor records remain separate from mature miRNAs.",
    "- Raw, RPM and name-map SHA-256 values match the HD-T1 output index.",
    "- Manifest SHA-256 matches the HD-T1 integrity report.",
    "",
    paste0("- ", integrity_dt$check, ": ", ifelse(integrity_dt$pass, "PASS", "FAIL"))
  )
)

run_end <- Sys.time()
pkg_names <- c("limma", "edgeR", "statmod", "data.table", "ggplot2", "digest")
pkg_versions <- data.table(
  package = pkg_names,
  version = vapply(pkg_names, function(p) as.character(packageVersion(p)), character(1))
)
write_csv(pkg_versions, file.path(result_root, "08_environment", "HD_T3_PACKAGE_VERSIONS.csv"))
capture.output(sessionInfo(), file = file.path(result_root, "08_environment", "HD_T3_SESSION_INFO.txt"))
write_md(
  file.path(result_root, "HD_T3_ENVIRONMENT_LOCK.md"),
  c(
    "# HD-T3 environment lock",
    "",
    paste0("- Start: ", format(run_start, "%Y-%m-%d %H:%M:%S %Z"), "."),
    paste0("- End: ", format(run_end, "%Y-%m-%d %H:%M:%S %Z"), "."),
    sprintf("- Runtime seconds: %.2f.", as.numeric(difftime(run_end, run_start, units = "secs"))),
    paste0("- R version: ", R.version.string, "."),
    paste0("- Operating system: ", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " | "), "."),
    paste0("- Project-local library: `", local_lib, "`."),
    "- Random seed: 20260716.",
    paste0("- Packages: ", paste(paste0(pkg_versions$package, " ", pkg_versions$version), collapse = "; "), "."),
    "- Plotting backend: R only.",
    "- HD-T1, HD-T2 and HD-R2 inputs remained read-only.",
    "- No external database or network resource was queried."
  )
)

write_md(
  file.path(result_root, "HD_T3_EXECUTIVE_SUMMARY.md"),
  c(
    "# HD-T3 executive summary",
    "",
    paste0("**Final state:** `", final_status, "`"),
    "",
    sprintf("- Primary-analysis patients: %d.", nrow(meta)),
    sprintf("- miRNAs before/after primary filtering: %d/%d.", nrow(raw_mat), sum(filter_counts$primary_after_n)),
    sprintf("- Mature miRNAs before/after: %d/%d.", filter_counts[feature_family == "MATURE", before_n], filter_counts[feature_family == "MATURE", primary_after_n]),
    sprintf("- Precursors before/after: %d/%d.", filter_counts[feature_family == "PRECURSOR", before_n], filter_counts[feature_family == "PRECURSOR", primary_after_n]),
    sprintf("- Adjusted mature-family FDR-significant miRNAs: %d.", primary_sig_mature),
    sprintf("- Unadjusted mature-family FDR-significant miRNAs: %d.", unadjusted_sig_mature),
    sprintf("- Reduced-zone mature-family FDR-significant miRNAs: %d.", reduced_sig_mature),
    sprintf("- Adjusted-stable mature miRNAs: %d.", stable_mature),
    sprintf("- Directionally stable low-precision mature miRNAs: %d.", low_precision_mature),
    sprintf("- Zone-dependent tested miRNAs: %d.", zone_dependent_n),
    sprintf("- Sample-sensitive tested miRNAs: %d.", sample_sensitive_n),
    sprintf("- Mature-family quality-weight range: %.4f to %.4f.", min(quality_weights$quality_weight), max(quality_weights$quality_weight)),
    paste0("- Independent miRNA information: `", independent_mirna_information, "`."),
    paste0("- HD-T4 eligibility: `", hd_t4_eligibility, "`."),
    "- HD-T4 automatically started: NO."
  )
)

write_md(
  file.path(result_root, "HD_T3_FINAL_GATE.md"),
  c(
    "# HD-T3 final gate",
    "",
    "## Final status",
    "",
    paste0("`", final_status, "`"),
    "",
    sprintf("- Final patients: %d.", nrow(meta)),
    sprintf("- miRNAs before/after filtering: %d/%d.", nrow(raw_mat), sum(filter_counts$primary_after_n)),
    sprintf("- Adjusted mature FDR-significant miRNAs: %d.", primary_sig_mature),
    sprintf("- Unadjusted mature FDR-significant miRNAs: %d.", unadjusted_sig_mature),
    sprintf("- Reduced-zone mature FDR-significant miRNAs: %d.", reduced_sig_mature),
    sprintf("- Adjusted-stable mature miRNAs: %d.", stable_mature),
    sprintf("- Low-precision stable mature miRNAs: %d.", low_precision_mature),
    sprintf("- Zone-dependent miRNAs: %d.", zone_dependent_n),
    sprintf("- Sample-sensitive miRNAs: %d.", sample_sensitive_n),
    sprintf("- Adjusted-stable precursor miRNAs: %d.", stable_precursor),
    sprintf("- Mature quality-weight range: %.6f to %.6f.", min(quality_weights$quality_weight), max(quality_weights$quality_weight)),
    paste0("- Tissue-zone effect: `", if (zone_dependent_n > 0L) "MATERIAL_FOR_A_SUBSET" else "NO_CLASSIFIED_ZONE_DEPENDENCE", "`."),
    paste0("- Independent of mRNA information: `", independent_mirna_information, "`."),
    paste0("- HD-T4 Work-review eligibility: `", hd_t4_eligibility, "`."),
    "- HD-T4, target analysis, miRNA-mRNA network, single-cell reference analysis and external validation were not started.",
    "- HARD STOP: wait for Work review."
  )
)
file.copy(
  file.path(result_root, "HD_T3_FINAL_GATE.md"),
  file.path(result_root, "09_gate", "HD_T3_FINAL_GATE.md"),
  overwrite = TRUE
)

log_msg("HD-T3 completed with status ", final_status, ". Hard stop enforced.")
file.copy(log_file, file.path(result_root, "HD_T3_ANALYSIS_LOG.md"), overwrite = TRUE)

all_outputs <- list.files(result_root, recursive = TRUE, full.names = TRUE)
all_outputs <- all_outputs[file.info(all_outputs)$isdir == FALSE]
hash_index <- data.table(
  relative_path = substring(
    normalizePath(all_outputs, winslash = "/", mustWork = TRUE),
    nchar(project_root) + 2L
  ),
  file_size_bytes = as.numeric(file.info(all_outputs)$size),
  sha256 = vapply(all_outputs, sha256, character(1))
)
hash_index <- hash_index[relative_path != "results/HD_T3/HD_T3_SHA256_INDEX.csv"]
write_csv(hash_index, file.path(result_root, "HD_T3_SHA256_INDEX.csv"))

cat(final_status, "\n")
