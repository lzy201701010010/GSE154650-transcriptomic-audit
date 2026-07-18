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

result_root <- file.path(project_root, "results", "HD_T2")
dirs <- c(
  "00_logs", "01_inputs", "02_gene_level", "03_transcript_level",
  "04_models", "05_sensitivity", "06_modules", "07_tables",
  "08_figures", "09_environment", "10_gate"
)
for (d in dirs) dir.create(file.path(result_root, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(result_root, "00_logs", "HD_T2_RUN.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
write_md <- function(path, lines) writeLines(enc2utf8(lines), path, useBytes = TRUE)
write_csv <- function(x, path) fwrite(x, path, quote = TRUE, na = "NA")
write_gz_csv <- function(x, path) fwrite(x, path, quote = TRUE, na = "NA", compress = "gzip")
write_gz_tsv <- function(x, path) fwrite(x, path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE)

run_start <- Sys.time()
log_msg("HD-T2 start; HD-T3 autostart prohibited.")

paths <- list(
  manifest = file.path(project_root, "HD_DATASET_AND_SAMPLE_MANIFEST.csv"),
  raw_enst = file.path(project_root, "results", "HD_T1", "02_matrices", "HD_T1_RAW_ENST.tsv.gz"),
  rpm_enst = file.path(project_root, "results", "HD_T1", "02_matrices", "HD_T1_RPM_ENST.tsv.gz"),
  annotation = file.path(project_root, "results", "HD_T1", "03_feature_annotation", "HD_T1_ENST_ANNOTATION_MAP.csv"),
  qc = file.path(project_root, "results", "HD_T1", "04_qc_tables", "HD_T1_SAMPLE_QC_METRICS.csv"),
  warnings = file.path(project_root, "results", "HD_T1", "04_qc_tables", "HD_T1_SAMPLE_QC_WARNINGS.csv"),
  t1_gate = file.path(project_root, "results", "HD_T1", "HD_T1_FINAL_GATE.md"),
  t1_input_manifest = file.path(project_root, "results", "HD_T1", "01_manifests", "HD_T1_INPUT_FILE_MANIFEST.csv"),
  module_registry = file.path(project_root, "source_metadata", "HD_T2_module_registry", "HD_T2_PREDEFINED_MODULE_REGISTRY_FROZEN.csv"),
  module_members = file.path(project_root, "source_metadata", "HD_T2_module_registry", "HD_T2_PREDEFINED_MODULE_MEMBERS_FROZEN.csv"),
  module_freeze_index = file.path(project_root, "source_metadata", "HD_T2_module_registry", "HD_T2_MODULE_FREEZE_INDEX.csv")
)
if (!all(file.exists(unlist(paths)))) {
  missing <- unlist(paths)[!file.exists(unlist(paths))]
  stop("Missing required inputs: ", paste(missing, collapse = "; "))
}

input_hashes <- data.table(
  input_name = names(paths),
  path = normalizePath(unlist(paths), winslash = "/", mustWork = TRUE),
  file_size_bytes = as.numeric(file.info(unlist(paths))$size),
  sha256 = vapply(unlist(paths), sha256, character(1))
)
write_csv(input_hashes, file.path(result_root, "01_inputs", "HD_T2_INPUT_SHA256.csv"))

manifest <- fread(paths$manifest)
raw_dt <- fread(paths$raw_enst)
rpm_dt <- fread(paths$rpm_enst)
annotation <- fread(paths$annotation)
qc <- fread(paths$qc)
warnings <- fread(paths$warnings)
module_registry <- fread(paths$module_registry)
module_members <- fread(paths$module_members)

required_manifest_cols <- c("patient_id", "GEO_sample_id", "disease_group", "sex", "BMI", "tissue_zone")
stopifnot(all(required_manifest_cols %in% names(manifest)))
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
integrity <- list(
  patients_38 = nrow(meta) == 38L,
  group_20_18 = identical(as.integer(table(meta$disease_group)), c(18L, 20L)),
  raw_dimensions = nrow(raw_dt) == 32153L && ncol(raw_dt) == 39L,
  rpm_dimensions = nrow(rpm_dt) == 32153L && ncol(rpm_dt) == 39L,
  columns_match_manifest = identical(patient_cols, meta$patient_id),
  rpm_columns_match_raw = identical(names(rpm_dt), names(raw_dt)),
  unique_features = !anyDuplicated(raw_dt$feature_id),
  annotation_unique_features = !anyDuplicated(annotation$original_feature_id),
  annotation_feature_set_matches_raw = setequal(annotation$original_feature_id, raw_dt$feature_id),
  no_duplicate_patients = !anyDuplicated(meta$patient_id),
  no_missing_covariates = !anyNA(meta[, .(disease_group, sex, BMI, tissue_zone)])
)

raw_mat <- as.matrix(raw_dt[, ..patient_cols])
rpm_mat <- as.matrix(rpm_dt[, ..patient_cols])
storage.mode(raw_mat) <- "double"
storage.mode(rpm_mat) <- "double"
rownames(raw_mat) <- raw_dt$feature_id
rownames(rpm_mat) <- rpm_dt$feature_id
integrity$raw_no_na <- !anyNA(raw_mat)
integrity$raw_no_inf <- all(is.finite(raw_mat))
integrity$raw_no_negative <- all(raw_mat >= 0)
integrity$rpm_no_na <- !anyNA(rpm_mat)
integrity$rpm_no_inf <- all(is.finite(rpm_mat))
integrity$rpm_no_negative <- all(rpm_mat >= 0)
integrity$fractional_precision_preserved <- mean(abs(raw_mat - round(raw_mat)) > 1e-8) > 0.5
integrity$warning_sample_not_zero_filled <- {
  j <- match("J01148", colnames(raw_mat))
  !is.na(j) && nrow(raw_mat) == 32153L && sum(raw_mat[, j] == 0) < nrow(raw_mat)
}

integrity_dt <- data.table(
  check = names(integrity),
  pass = vapply(integrity, isTRUE, logical(1))
)
write_csv(integrity_dt, file.path(result_root, "01_inputs", "HD_T2_INPUT_INTEGRITY_CHECKS.csv"))
if (!all(integrity_dt$pass)) {
  stop("HD-T2 input integrity failure: ", paste(integrity_dt[pass == FALSE, check], collapse = ", "))
}
annotation <- annotation[match(raw_dt$feature_id, original_feature_id)]
stopifnot(identical(annotation$original_feature_id, raw_dt$feature_id))
log_msg("Input integrity passed: 32153 ENST x 38 patients; 20 cases and 18 controls.")

design_primary <- model.matrix(~ disease_group + sex + BMI_centered + tissue_zone, data = meta)
design_unadjusted <- model.matrix(~ disease_group, data = meta)
design_no_zone <- model.matrix(~ disease_group + sex + BMI_centered, data = meta)
if (qr(design_primary)$rank != ncol(design_primary)) stop("Primary design is not full rank.")

mapped <- annotation$mapped_gene_n == 1L &
  !is.na(annotation$ensembl_gene_id) &
  nzchar(annotation$ensembl_gene_id)
annotation[, primary_gene_inclusion := ifelse(mapped, "INCLUDED_UNIQUE_GENE_MAPPING", "EXCLUDED_UNMAPPED_OR_RETIRED")]
annotation[, transcript_n_per_gene := fifelse(
  mapped,
  ave(as.integer(mapped), ensembl_gene_id, FUN = sum),
  NA_integer_
)]
annotation[, duplicate_gene_symbol_across_ids := FALSE]
write_csv(annotation, file.path(result_root, "HD_T2_GENE_AGGREGATION_MAP.csv"))
write_csv(
  annotation[!mapped],
  file.path(result_root, "02_gene_level", "HD_T2_UNMAPPED_AND_AMBIGUOUS_ENST.csv")
)

gene_raw <- rowsum(raw_mat[mapped, , drop = FALSE], annotation$ensembl_gene_id[mapped], reorder = FALSE)
gene_rpm <- rowsum(rpm_mat[mapped, , drop = FALSE], annotation$ensembl_gene_id[mapped], reorder = FALSE)
gene_meta <- unique(annotation[mapped, .(
  ensembl_gene_id, gene_symbol, gene_biotype, mapping_status
)], by = "ensembl_gene_id")
gene_meta <- gene_meta[match(rownames(gene_raw), ensembl_gene_id)]
stopifnot(identical(gene_meta$ensembl_gene_id, rownames(gene_raw)))

gene_raw_out <- data.table(ensembl_gene_id = rownames(gene_raw), as.data.table(gene_raw))
gene_rpm_out <- data.table(ensembl_gene_id = rownames(gene_rpm), as.data.table(gene_rpm))
write_gz_tsv(gene_raw_out, file.path(result_root, "02_gene_level", "HD_T2_GENE_LEVEL_RAW_MATRIX.tsv.gz"))
write_gz_tsv(gene_rpm_out, file.path(result_root, "02_gene_level", "HD_T2_GENE_LEVEL_RPM_DESCRIPTIVE.tsv.gz"))

aggregation_summary <- c(
  "# HD-T2 gene aggregation summary",
  "",
  sprintf("- Input ENST: %d.", nrow(annotation)),
  sprintf("- Uniquely mapped ENST: %d.", sum(mapped)),
  sprintf("- ENST entering primary gene aggregation: %d.", sum(mapped)),
  sprintf("- Final Ensembl genes: %d.", nrow(gene_raw)),
  sprintf("- Unmapped or retired ENST: %d.", sum(!mapped)),
  sprintf("- One-to-many ENST mappings: %d.", sum(annotation$mapped_gene_n > 1L, na.rm = TRUE)),
  sprintf("- Multi-transcript genes: %d.", sum(table(annotation$ensembl_gene_id[mapped]) > 1L)),
  "- Aggregation: sum fractional raw abundance across uniquely mapped ENST sharing one Ensembl gene ID.",
  "- Gene symbols are annotations only and were not used as aggregation keys.",
  "- No absent feature was filled with zero and no fractional raw value was rounded."
)
write_md(file.path(result_root, "HD_T2_GENE_AGGREGATION_SUMMARY.md"), aggregation_summary)
log_msg("Gene aggregation completed: ", nrow(gene_raw), " genes from ", sum(mapped), " uniquely mapped ENST.")

dge_gene_all <- DGEList(gene_raw)
gene_keep_primary <- filterByExpr(dge_gene_all, design = design_primary)
gene_keep_relaxed <- filterByExpr(
  dge_gene_all, design = design_primary, min.count = 5, min.total.count = 10
)
gene_keep_stringent <- filterByExpr(
  dge_gene_all, design = design_primary, min.count = 20, min.total.count = 30
)
dge_tx_all <- DGEList(raw_mat)
tx_keep_primary <- filterByExpr(dge_tx_all, design = design_primary)

gene_filter_audit <- copy(gene_meta)
gene_filter_audit[, `:=`(
  total_raw = rowSums(gene_raw),
  nonzero_patient_n = rowSums(gene_raw > 0),
  keep_primary_filterByExpr = gene_keep_primary,
  keep_relaxed_min5_total10 = gene_keep_relaxed,
  keep_stringent_min20_total30 = gene_keep_stringent
)]
tx_filter_audit <- annotation[, .(
  original_feature_id, source_enst_versioned, ensembl_gene_id, gene_symbol,
  mapped_gene_n, mapping_status, primary_gene_inclusion
)]
tx_filter_audit[, `:=`(
  total_raw = rowSums(raw_mat),
  nonzero_patient_n = rowSums(raw_mat > 0),
  keep_primary_filterByExpr = tx_keep_primary
)]
write_csv(gene_filter_audit, file.path(result_root, "HD_T2_GENE_FILTER_AUDIT.csv"))
write_csv(tx_filter_audit, file.path(result_root, "HD_T2_TRANSCRIPT_FILTER_AUDIT.csv"))

filter_summary <- c(
  "# HD-T2 filtering summary",
  "",
  sprintf("- Genes before filtering: %d.", length(gene_keep_primary)),
  sprintf("- Genes after primary filterByExpr: %d.", sum(gene_keep_primary)),
  sprintf("- Genes after relaxed sensitivity filter: %d.", sum(gene_keep_relaxed)),
  sprintf("- Genes after stringent sensitivity filter: %d.", sum(gene_keep_stringent)),
  sprintf("- ENST before filtering: %d.", length(tx_keep_primary)),
  sprintf("- ENST after primary filterByExpr: %d.", sum(tx_keep_primary)),
  "- Primary filtering used the frozen adjusted design and package defaults.",
  "- Relaxed and stringent thresholds were fixed in HD_T2_PROMPT_REVIEW.md before outcome inspection.",
  "- No feature was manually retained because of biological interest."
)
write_md(file.path(result_root, "HD_T2_FILTERING_SUMMARY.md"), filter_summary)
log_msg("Filtering completed: genes ", sum(gene_keep_primary), "/", length(gene_keep_primary),
        "; ENST ", sum(tx_keep_primary), "/", length(tx_keep_primary), ".")

gene_counts <- gene_raw[gene_keep_primary, , drop = FALSE]
gene_meta_tested <- gene_meta[gene_keep_primary]
tx_counts <- raw_mat[tx_keep_primary, , drop = FALSE]
tx_annotation_tested <- annotation[tx_keep_primary]

coef_name <- "disease_grouphemorrhoidal disease"
fit_voom_model <- function(counts, model_meta, formula, model_name, save_plot = FALSE) {
  design <- model.matrix(formula, data = model_meta)
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    return(list(
      status = "NOT_ESTIMABLE", model_name = model_name, design = design,
      rank = rank, result = NULL, v = NULL, fit = NULL, coef_index = NA_integer_
    ))
  }
  coef_index <- match(coef_name, colnames(design))
  if (is.na(coef_index)) {
    return(list(
      status = "DISEASE_COEFFICIENT_NOT_FOUND", model_name = model_name,
      design = design, rank = rank, result = NULL, v = NULL, fit = NULL,
      coef_index = NA_integer_
    ))
  }
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMM")
  v <- voomWithQualityWeights(
    y, design = design, plot = FALSE, save.plot = save_plot
  )
  fit <- eBayes(lmFit(v, design), robust = FALSE)
  tt <- topTable(fit, coef = coef_index, number = Inf, sort.by = "none", confint = 0.95)
  result <- data.table(
    feature_id = rownames(tt),
    log2_fold_change = tt$logFC,
    average_expression = tt$AveExpr,
    moderated_t = tt$t,
    p_value = tt$P.Value,
    BH_FDR = tt$adj.P.Val,
    B_statistic = tt$B,
    CI_95_low = tt$CI.L,
    CI_95_high = tt$CI.R
  )
  result[, moderated_standard_error := fifelse(
    is.finite(moderated_t) & abs(moderated_t) > 0,
    abs(log2_fold_change / moderated_t),
    NA_real_
  )]
  sw <- v$targets$sample.weights
  if (is.null(sw)) {
    sw <- apply(v$weights, 2, median, na.rm = TRUE)
    sw <- sw / exp(mean(log(sw)))
  }
  names(sw) <- colnames(counts)
  kish_n <- (sum(sw)^2) / sum(sw^2)
  result[, `:=`(
    n_patients = nrow(model_meta),
    effective_sample_size_kish = kish_n,
    model_name = model_name,
    statistical_unit = "patient/library",
    contrast = "hemorrhoidal disease minus anal fissure surgical control"
  )]
  list(
    status = "ESTIMABLE", model_name = model_name, design = design, rank = rank,
    result = result, v = v, fit = fit, coef_index = coef_index,
    sample_weights = sw, norm_factors = y$samples$norm.factors
  )
}

attach_gene_annotation <- function(res) {
  out <- merge(
    res,
    gene_meta_tested,
    by.x = "feature_id",
    by.y = "ensembl_gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  setnames(out, "feature_id", "ensembl_gene_id")
  setcolorder(out, c(
    "ensembl_gene_id", "gene_symbol", "gene_biotype", "mapping_status",
    setdiff(names(out), c("ensembl_gene_id", "gene_symbol", "gene_biotype", "mapping_status"))
  ))
  out
}

model_primary <- fit_voom_model(
  gene_counts, meta,
  ~ disease_group + sex + BMI_centered + tissue_zone,
  "PRIMARY_ADJUSTED", save_plot = TRUE
)
model_unadjusted <- fit_voom_model(
  gene_counts, meta, ~ disease_group, "UNADJUSTED"
)
model_no_zone <- fit_voom_model(
  gene_counts, meta,
  ~ disease_group + sex + BMI_centered,
  "NO_ZONE_ADJUSTMENT"
)
if (model_primary$status != "ESTIMABLE") stop("Primary model failed.")

reduced_idx <- meta$tissue_zone != "anoderm region"
model_reduced <- fit_voom_model(
  gene_counts[, reduced_idx, drop = FALSE],
  droplevels(meta[reduced_idx]),
  ~ disease_group + sex + BMI_centered + tissue_zone,
  "REDUCED_ZONE_EXCLUDING_ANODERM"
)
transition_idx <- meta$tissue_zone == "transition zone"
model_transition <- fit_voom_model(
  gene_counts[, transition_idx, drop = FALSE],
  droplevels(meta[transition_idx]),
  ~ disease_group + sex + BMI_centered,
  "TRANSITION_ZONE_LOW_PRECISION"
)
intestinal_idx <- meta$tissue_zone == "intestinal mucosa region"
model_intestinal <- fit_voom_model(
  gene_counts[, intestinal_idx, drop = FALSE],
  droplevels(meta[intestinal_idx]),
  ~ disease_group + sex + BMI_centered,
  "INTESTINAL_MUCOSA_ZONE_LOW_PRECISION"
)

model_objects <- list(
  primary = model_primary,
  unadjusted = model_unadjusted,
  no_zone = model_no_zone,
  reduced = model_reduced,
  transition = model_transition,
  intestinal = model_intestinal
)
design_audit <- rbindlist(lapply(model_objects, function(x) data.table(
  model_name = x$model_name,
  status = x$status,
  n_patients = nrow(x$design),
  design_columns = ncol(x$design),
  design_rank = x$rank,
  residual_degrees_of_freedom = nrow(x$design) - x$rank,
  coefficient_names = paste(colnames(x$design), collapse = ";")
)), fill = TRUE)
write_csv(design_audit, file.path(result_root, "04_models", "HD_T2_DESIGN_MATRIX_AUDIT.csv"))

formal_gene_outputs <- list(
  HD_T2_GENE_DE_ADJUSTED.csv.gz = model_primary,
  HD_T2_GENE_DE_UNADJUSTED.csv.gz = model_unadjusted,
  HD_T2_GENE_DE_NO_ZONE_ADJUSTMENT.csv.gz = model_no_zone,
  HD_T2_GENE_DE_REDUCED_ZONE.csv.gz = model_reduced,
  HD_T2_GENE_DE_TRANSITION_ZONE.csv.gz = model_transition,
  HD_T2_GENE_DE_INTESTINAL_MUCOSA_ZONE.csv.gz = model_intestinal
)
for (nm in names(formal_gene_outputs)) {
  obj <- formal_gene_outputs[[nm]]
  if (obj$status == "ESTIMABLE") {
    write_gz_csv(attach_gene_annotation(obj$result), file.path(result_root, nm))
  } else {
    write_gz_csv(data.table(
      model_name = obj$model_name, status = obj$status,
      statistical_unit = "patient/library"
    ), file.path(result_root, nm))
  }
}
log_msg("Primary and mandatory zone models completed.")

model_spec <- c(
  "# HD-T2 model specification",
  "",
  "- Primary estimand: adjusted patient-level log2 expression difference for hemorrhoidal disease minus anal fissure surgical control.",
  "- Primary model: `~ disease_group + sex + BMI_centered + tissue_zone`.",
  "- Reference levels: anal fissure surgical control; female; anoderm region.",
  "- Input: fractional raw abundance aggregated by Ensembl gene ID without rounding.",
  "- Filtering: edgeR::filterByExpr with the primary adjusted design.",
  "- Normalization: TMM.",
  "- Mean-variance model: limma::voomWithQualityWeights.",
  "- Inference: limma::lmFit followed by limma::eBayes(robust = FALSE).",
  "- Multiplicity: Benjamini-Hochberg within each formal result family.",
  "- Tissue-zone interaction: not run because higher-priority frozen plans did not explicitly authorize it.",
  "- Age and batch: excluded because structurally unavailable.",
  "- All formal case-control models use the patient/library as the biological statistical unit.",
  "",
  "## Stability rules fixed before outcome inspection",
  "",
  "- Near-zero classification threshold: |log2FC| < 0.05.",
  "- Sample-sensitive: leave-one-out direction retention <80%, warning-sample exclusion direction reversal, or leave-one-out maximum absolute change greater than max(0.5, |primary log2FC|).",
  "- Zone-dependent: direction reversal after zone removal/restriction or at least 50% attenuation when |primary log2FC| >=0.10.",
  "- Robust: primary FDR <0.05, consistent non-zero direction across primary, unadjusted, no-zone, reduced-zone and warning-sample exclusion analyses, with leave-one-out retention >=90%.",
  "- Directionally stable low precision: consistent direction with leave-one-out retention >=80% but not meeting the robust FDR rule.",
  "- Transcript heterogeneity: at least two tested transcripts with |log2FC| >=0.10 in opposite directions."
)
write_md(file.path(result_root, "HD_T2_MODEL_SPECIFICATION.md"), model_spec)

warning_exclusions <- list(
  EXCLUDE_GSM4676447 = "J01138",
  EXCLUDE_GSM4676457 = "J01148",
  EXCLUDE_BOTH_WARNING_FILES = c("J01138", "J01148")
)
warning_models <- list()
for (nm in names(warning_exclusions)) {
  keep_patients <- !meta$patient_id %chin% warning_exclusions[[nm]]
  warning_models[[nm]] <- fit_voom_model(
    gene_counts[, keep_patients, drop = FALSE],
    droplevels(meta[keep_patients]),
    ~ disease_group + sex + BMI_centered + tissue_zone,
    nm
  )
  if (warning_models[[nm]]$status == "ESTIMABLE") {
    write_gz_csv(
      attach_gene_annotation(warning_models[[nm]]$result),
      file.path(result_root, "05_sensitivity", paste0("HD_T2_GENE_DE_", nm, ".csv.gz"))
    )
  }
}

fit_filter_sensitivity <- function(keep, name) {
  counts <- gene_raw[keep, , drop = FALSE]
  gm <- gene_meta[keep]
  obj <- fit_voom_model(
    counts, meta,
    ~ disease_group + sex + BMI_centered + tissue_zone,
    name
  )
  if (obj$status == "ESTIMABLE") {
    out <- merge(
      obj$result, gm,
      by.x = "feature_id", by.y = "ensembl_gene_id",
      all.x = TRUE, sort = FALSE
    )
    write_gz_csv(out, file.path(result_root, "05_sensitivity", paste0("HD_T2_", name, ".csv.gz")))
  }
  obj
}
filter_relaxed_model <- fit_filter_sensitivity(gene_keep_relaxed, "FILTER_RELAXED_MIN5_TOTAL10")
filter_stringent_model <- fit_filter_sensitivity(gene_keep_stringent, "FILTER_STRINGENT_MIN20_TOTAL30")

sex_models <- list()
for (sex_level in levels(meta$sex)) {
  idx <- meta$sex == sex_level
  sex_models[[sex_level]] <- fit_voom_model(
    gene_counts[, idx, drop = FALSE],
    droplevels(meta[idx]),
    ~ disease_group + BMI_centered + tissue_zone,
    paste0("SEX_STRATIFIED_", toupper(sex_level))
  )
  if (sex_models[[sex_level]]$status == "ESTIMABLE") {
    write_gz_csv(
      attach_gene_annotation(sex_models[[sex_level]]$result),
      file.path(result_root, "05_sensitivity", paste0("HD_T2_GENE_DE_SEX_", toupper(sex_level), ".csv.gz"))
    )
  }
}

fit_interaction <- function(counts, model_meta) {
  design <- model.matrix(~ disease_group * sex + BMI_centered + tissue_zone, data = model_meta)
  interaction_name <- grep("disease_group.*:sex", colnames(design), value = TRUE)
  if (qr(design)$rank != ncol(design) || length(interaction_name) != 1L) {
    return(list(status = "NOT_ESTIMABLE", design = design, result = NULL))
  }
  y <- calcNormFactors(DGEList(counts), method = "TMM")
  v <- voomWithQualityWeights(y, design = design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = FALSE)
  tt <- topTable(fit, coef = interaction_name, number = Inf, sort.by = "none", confint = 0.95)
  list(
    status = "ESTIMABLE",
    design = design,
    result = data.table(
      ensembl_gene_id = rownames(tt),
      interaction_log2_effect = tt$logFC,
      average_expression = tt$AveExpr,
      moderated_t = tt$t,
      p_value = tt$P.Value,
      BH_FDR = tt$adj.P.Val,
      CI_95_low = tt$CI.L,
      CI_95_high = tt$CI.R,
      model_name = "EXPLORATORY_DISEASE_BY_SEX_INTERACTION",
      statistical_unit = "patient/library"
    )
  )
}
sex_interaction <- fit_interaction(gene_counts, meta)
if (sex_interaction$status == "ESTIMABLE") {
  sex_interaction$result <- merge(sex_interaction$result, gene_meta_tested, by = "ensembl_gene_id", all.x = TRUE)
  write_gz_csv(
    sex_interaction$result,
    file.path(result_root, "05_sensitivity", "HD_T2_DISEASE_BY_SEX_INTERACTION.csv.gz")
  )
}

fractional_fraction <- mean(abs(raw_mat - round(raw_mat)) > 1e-8)
edger_eligibility <- data.table(
  method = "edgeR quasi-likelihood robust dispersion",
  fractional_raw_value_fraction = fractional_fraction,
  rounding_or_coercion_allowed = "NO",
  status = if (fractional_fraction > 0.5) "NOT_RUN_FRACTIONAL_INPUT_INELIGIBLE" else "ELIGIBLE",
  rationale = if (fractional_fraction > 0.5)
    "Deposited raw abundance is predominantly fractional; count-likelihood sensitivity was not forced."
  else
    "Input diagnostic passed."
)
write_csv(edger_eligibility, file.path(result_root, "05_sensitivity", "HD_T2_EDGER_QL_ELIGIBILITY.csv"))

primary_res <- copy(model_primary$result)
setnames(primary_res, "feature_id", "ensembl_gene_id")
unadj_res <- copy(model_unadjusted$result)
setnames(unadj_res, "feature_id", "ensembl_gene_id")
nozone_res <- copy(model_no_zone$result)
setnames(nozone_res, "feature_id", "ensembl_gene_id")
reduced_res <- if (model_reduced$status == "ESTIMABLE") {
  x <- copy(model_reduced$result)
  setnames(x, "feature_id", "ensembl_gene_id")
  x
} else NULL
transition_res <- if (model_transition$status == "ESTIMABLE") {
  x <- copy(model_transition$result)
  setnames(x, "feature_id", "ensembl_gene_id")
  x
} else NULL
intestinal_res <- if (model_intestinal$status == "ESTIMABLE") {
  x <- copy(model_intestinal$result)
  setnames(x, "feature_id", "ensembl_gene_id")
  x
} else NULL

warn_gene <- primary_res[, .(
  ensembl_gene_id,
  primary_log2FC = log2_fold_change,
  primary_FDR = BH_FDR
)]
scenario_summary <- list()
for (nm in names(warning_models)) {
  obj <- warning_models[[nm]]
  if (obj$status != "ESTIMABLE") next
  tmp <- copy(obj$result)
  setnames(tmp, c("feature_id", "log2_fold_change", "BH_FDR"),
           c("ensembl_gene_id", paste0(nm, "_log2FC"), paste0(nm, "_FDR")))
  tmp <- tmp[, c("ensembl_gene_id", paste0(nm, "_log2FC"), paste0(nm, "_FDR")), with = FALSE]
  warn_gene <- merge(warn_gene, tmp, by = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
  scenario_summary[[nm]] <- data.table(
    scenario = nm,
    n_patients = nrow(obj$design),
    significant_gene_n = sum(obj$result$BH_FDR < 0.05, na.rm = TRUE),
    log2FC_spearman_vs_primary = cor(
      obj$result$log2_fold_change, model_primary$result$log2_fold_change,
      method = "spearman", use = "complete.obs"
    ),
    direction_concordance_vs_primary = mean(
      sign(obj$result$log2_fold_change) == sign(model_primary$result$log2_fold_change),
      na.rm = TRUE
    ),
    median_absolute_log2FC_change = median(
      abs(obj$result$log2_fold_change - model_primary$result$log2_fold_change),
      na.rm = TRUE
    )
  )
}
warn_gene <- merge(warn_gene, gene_meta_tested, by = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
write_csv(warn_gene, file.path(result_root, "HD_T2_WARNING_SAMPLE_SENSITIVITY.csv"))
write_csv(rbindlist(scenario_summary), file.path(result_root, "05_sensitivity", "HD_T2_WARNING_SAMPLE_SCENARIO_SUMMARY.csv"))

module_index_from_v <- function(v) {
  ids <- rownames(v$E)
  symbols <- gene_meta_tested$gene_symbol[match(ids, gene_meta_tested$ensembl_gene_id)]
  setNames(lapply(module_registry$module_id, function(mid) {
    genes <- module_members[module_id == mid, unique(source_gene_symbol)]
    which(!is.na(symbols) & symbols %chin% genes)
  }), module_registry$module_id)
}
camera_for_model <- function(obj) {
  if (obj$status != "ESTIMABLE") return(NULL)
  idx <- module_index_from_v(obj$v)
  cam <- camera(
    obj$v, index = idx, design = obj$design,
    contrast = obj$coef_index, inter.gene.cor = 0.01, sort = FALSE
  )
  out <- data.table(module_id = rownames(cam), cam)
  setnames(out, c("NGenes", "Direction", "PValue", "FDR"),
           c("tested_gene_n", "camera_direction", "camera_p_value", "camera_BH_FDR"))
  out
}

module_camera <- list(
  PRIMARY_ADJUSTED = camera_for_model(model_primary),
  UNADJUSTED = camera_for_model(model_unadjusted),
  NO_ZONE_ADJUSTMENT = camera_for_model(model_no_zone),
  REDUCED_ZONE = camera_for_model(model_reduced),
  EXCLUDE_BOTH_WARNING_FILES = camera_for_model(warning_models$EXCLUDE_BOTH_WARNING_FILES)
)

primary_idx <- module_index_from_v(model_primary$v)
z_expression <- t(scale(t(model_primary$v$E)))
z_expression[!is.finite(z_expression)] <- NA_real_
module_scores <- sapply(primary_idx, function(idx) {
  if (length(idx) == 0L) return(rep(NA_real_, ncol(z_expression)))
  colMeans(z_expression[idx, , drop = FALSE], na.rm = TRUE)
})
module_scores <- t(module_scores)
rownames(module_scores) <- names(primary_idx)
colnames(module_scores) <- colnames(model_primary$v$E)
module_score_out <- data.table(
  module_id = rep(rownames(module_scores), times = ncol(module_scores)),
  patient_id = rep(colnames(module_scores), each = nrow(module_scores)),
  module_z_score = as.vector(module_scores)
)
module_score_out <- merge(module_score_out, meta, by = "patient_id", all.x = TRUE, sort = FALSE)
write_csv(module_score_out, file.path(result_root, "HD_T2_MODULE_SCORE_MATRIX.csv"))

module_score_fit <- eBayes(lmFit(module_scores, model_primary$design), robust = FALSE)
module_score_tt <- topTable(
  module_score_fit, coef = model_primary$coef_index,
  number = Inf, sort.by = "none", confint = 0.95
)
module_score_res <- data.table(
  module_id = rownames(module_score_tt),
  score_log2_difference = module_score_tt$logFC,
  score_moderated_t = module_score_tt$t,
  score_p_value = module_score_tt$P.Value,
  score_BH_FDR = module_score_tt$adj.P.Val,
  score_CI_95_low = module_score_tt$CI.L,
  score_CI_95_high = module_score_tt$CI.R
)

module_manifest <- copy(module_registry)
mapped_counts <- data.table(
  module_id = names(primary_idx),
  mapped_tested_gene_n = lengths(primary_idx)
)
module_manifest <- merge(module_manifest, mapped_counts, by = "module_id", all.x = TRUE, sort = FALSE)
module_manifest[, mapping_success_fraction := mapped_tested_gene_n / source_gene_n]
write_csv(module_manifest, file.path(result_root, "HD_T2_PREDEFINED_MODULE_MANIFEST.csv"))

module_diff <- merge(module_camera$PRIMARY_ADJUSTED, module_score_res, by = "module_id", all = TRUE)
module_diff <- merge(
  module_diff,
  module_manifest[, .(module_id, module_name, module_role, source_library, source_term, source_version)],
  by = "module_id", all.x = TRUE, sort = FALSE
)
module_diff[, `:=`(
  model_name = "PRIMARY_ADJUSTED",
  statistical_unit = "patient/library"
)]
write_csv(module_diff, file.path(result_root, "HD_T2_MODULE_DIFFERENTIAL_RESULTS.csv"))

loo_logfc <- matrix(
  NA_real_, nrow = nrow(gene_counts), ncol = nrow(meta),
  dimnames = list(rownames(gene_counts), meta$patient_id)
)
loo_fdr <- matrix(
  NA_real_, nrow = nrow(gene_counts), ncol = nrow(meta),
  dimnames = list(rownames(gene_counts), meta$patient_id)
)
loo_module <- list()
sample_influence <- vector("list", nrow(meta))
log_msg("Starting leave-one-patient-out refits: ", nrow(meta), " iterations.")
for (i in seq_len(nrow(meta))) {
  keep <- seq_len(nrow(meta)) != i
  obj <- fit_voom_model(
    gene_counts[, keep, drop = FALSE],
    droplevels(meta[keep]),
    ~ disease_group + sex + BMI_centered + tissue_zone,
    paste0("LOO_", meta$patient_id[i])
  )
  if (obj$status == "ESTIMABLE") {
    loo_logfc[, i] <- obj$result$log2_fold_change
    loo_fdr[, i] <- obj$result$BH_FDR
    cam <- camera_for_model(obj)
    if (!is.null(cam)) {
      cam[, omitted_patient_id := meta$patient_id[i]]
      loo_module[[i]] <- cam
    }
    delta <- obj$result$log2_fold_change - model_primary$result$log2_fold_change
    sample_influence[[i]] <- data.table(
      patient_id = meta$patient_id[i],
      GEO_sample_id = meta$GEO_sample_id[i],
      model_status = obj$status,
      median_absolute_log2FC_change = median(abs(delta), na.rm = TRUE),
      q95_absolute_log2FC_change = quantile(abs(delta), 0.95, na.rm = TRUE),
      maximum_absolute_log2FC_change = max(abs(delta), na.rm = TRUE),
      direction_flip_fraction = mean(
        sign(obj$result$log2_fold_change) != sign(model_primary$result$log2_fold_change),
        na.rm = TRUE
      ),
      rank_spearman_vs_primary = cor(
        rank(abs(obj$result$moderated_t), ties.method = "average"),
        rank(abs(model_primary$result$moderated_t), ties.method = "average"),
        method = "spearman", use = "complete.obs"
      )
    )
  } else {
    sample_influence[[i]] <- data.table(
      patient_id = meta$patient_id[i],
      GEO_sample_id = meta$GEO_sample_id[i],
      model_status = obj$status
    )
  }
  if (i %% 5L == 0L || i == nrow(meta)) {
    log_msg("LOO progress ", i, "/", nrow(meta), ".")
  }
}

near_zero <- function(x) abs(x) < 0.05
primary_vec <- model_primary$result$log2_fold_change
loo_same <- sweep(loo_logfc, 1, primary_vec, function(x, p) {
  (near_zero(x) & near_zero(p)) | (!near_zero(x) & !near_zero(p) & sign(x) == sign(p))
})
loo_summary <- data.table(
  ensembl_gene_id = rownames(loo_logfc),
  primary_log2FC = primary_vec,
  direction_retention_fraction = rowMeans(loo_same, na.rm = TRUE),
  loo_min_log2FC = apply(loo_logfc, 1, min, na.rm = TRUE),
  loo_max_log2FC = apply(loo_logfc, 1, max, na.rm = TRUE),
  loo_max_absolute_change = apply(abs(loo_logfc - primary_vec), 1, max, na.rm = TRUE),
  loo_median_absolute_change = apply(abs(loo_logfc - primary_vec), 1, median, na.rm = TRUE),
  loo_FDR_lt_0_05_fraction = rowMeans(loo_fdr < 0.05, na.rm = TRUE)
)
loo_summary <- merge(loo_summary, gene_meta_tested, by = "ensembl_gene_id", all.x = TRUE, sort = FALSE)
write_csv(loo_summary, file.path(result_root, "HD_T2_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv"))
write_gz_csv(
  data.table(ensembl_gene_id = rownames(loo_logfc), as.data.table(loo_logfc)),
  file.path(result_root, "05_sensitivity", "HD_T2_LEAVE_ONE_PATIENT_OUT_LOG2FC_MATRIX.csv.gz")
)

sample_influence_dt <- rbindlist(sample_influence, fill = TRUE)
quality_weights <- data.table(
  patient_id = names(model_primary$sample_weights),
  quality_weight = as.numeric(model_primary$sample_weights),
  normalization_factor = as.numeric(model_primary$norm_factors)
)
quality_weights <- merge(quality_weights, meta, by = "patient_id", all.x = TRUE, sort = FALSE)
quality_weights <- merge(
  quality_weights,
  warnings[, .(patient_id, warning_class, warning_reasons, primary_analysis_status)],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
quality_weights <- merge(
  quality_weights,
  qc[analyte == "ENST", .(
    patient_id, median_spearman_to_others, minimum_spearman_to_others,
    PCA_distance, nonzero_feature_n, detection_fraction
  )],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
write_csv(quality_weights, file.path(result_root, "HD_T2_SAMPLE_QUALITY_WEIGHTS.csv"))

sample_influence_dt <- merge(sample_influence_dt, quality_weights, by = c("patient_id", "GEO_sample_id"), all.x = TRUE)
sample_influence_dt[, leverage_primary_design := hat(model_primary$design)]
write_csv(sample_influence_dt, file.path(result_root, "HD_T2_SAMPLE_INFLUENCE_METRICS.csv"))

weight_association <- list(
  disease = summary(lm(log(quality_weight) ~ disease_group, data = quality_weights)),
  zone = summary(lm(log(quality_weight) ~ tissue_zone, data = quality_weights))
)
capture.output(
  weight_association,
  file = file.path(result_root, "05_sensitivity", "HD_T2_QUALITY_WEIGHT_ASSOCIATION_MODELS.txt")
)

comparison <- Reduce(function(x, y) merge(x, y, by = "ensembl_gene_id", all = TRUE, sort = FALSE), list(
  primary_res[, .(ensembl_gene_id, primary_logFC = log2_fold_change, primary_FDR = BH_FDR, primary_t = moderated_t)],
  unadj_res[, .(ensembl_gene_id, unadjusted_logFC = log2_fold_change, unadjusted_FDR = BH_FDR, unadjusted_t = moderated_t)],
  nozone_res[, .(ensembl_gene_id, no_zone_logFC = log2_fold_change, no_zone_FDR = BH_FDR, no_zone_t = moderated_t)],
  if (!is.null(reduced_res)) reduced_res[, .(ensembl_gene_id, reduced_zone_logFC = log2_fold_change, reduced_zone_FDR = BH_FDR)] else primary_res[, .(ensembl_gene_id, reduced_zone_logFC = NA_real_, reduced_zone_FDR = NA_real_)],
  if (!is.null(transition_res)) transition_res[, .(ensembl_gene_id, transition_zone_logFC = log2_fold_change, transition_zone_FDR = BH_FDR)] else primary_res[, .(ensembl_gene_id, transition_zone_logFC = NA_real_, transition_zone_FDR = NA_real_)],
  if (!is.null(intestinal_res)) intestinal_res[, .(ensembl_gene_id, intestinal_zone_logFC = log2_fold_change, intestinal_zone_FDR = BH_FDR)] else primary_res[, .(ensembl_gene_id, intestinal_zone_logFC = NA_real_, intestinal_zone_FDR = NA_real_)],
  warn_gene[, .(
    ensembl_gene_id,
    exclude_both_logFC = EXCLUDE_BOTH_WARNING_FILES_log2FC,
    exclude_both_FDR = EXCLUDE_BOTH_WARNING_FILES_FDR
  )],
  loo_summary[, .(
    ensembl_gene_id, direction_retention_fraction,
    loo_max_absolute_change, loo_median_absolute_change,
    loo_FDR_lt_0_05_fraction
  )]
))
comparison <- merge(comparison, gene_meta_tested, by = "ensembl_gene_id", all.x = TRUE, sort = FALSE)

direction_agree <- function(a, b) {
  (near_zero(a) & near_zero(b)) | (!near_zero(a) & !near_zero(b) & sign(a) == sign(b))
}
comparison[, primary_rank := frank(-abs(primary_t), ties.method = "average")]
comparison[, unadjusted_rank := frank(-abs(unadjusted_t), ties.method = "average")]
comparison[, no_zone_rank := frank(-abs(no_zone_t), ties.method = "average")]
comparison[, rank_range := pmax(primary_rank, unadjusted_rank, no_zone_rank, na.rm = TRUE) -
  pmin(primary_rank, unadjusted_rank, no_zone_rank, na.rm = TRUE)]
comparison[, FDR_significant_model_n := rowSums(cbind(
  primary_FDR < 0.05, unadjusted_FDR < 0.05, no_zone_FDR < 0.05,
  reduced_zone_FDR < 0.05, exclude_both_FDR < 0.05
), na.rm = TRUE)]

comparison[, sample_sensitive := (
  direction_retention_fraction < 0.80 |
    !direction_agree(primary_logFC, exclude_both_logFC) |
    loo_max_absolute_change > pmax(0.5, abs(primary_logFC))
)]
comparison[, zone_dependent := (
  !direction_agree(primary_logFC, no_zone_logFC) |
    !direction_agree(primary_logFC, reduced_zone_logFC) |
    (
      abs(primary_logFC) >= 0.10 &
        abs(reduced_zone_logFC) < 0.50 * abs(primary_logFC)
    )
)]
comparison[, model_dependent := (
  !direction_agree(primary_logFC, unadjusted_logFC) |
    !direction_agree(primary_logFC, no_zone_logFC)
)]
comparison[, all_core_direction_stable := (
  direction_agree(primary_logFC, unadjusted_logFC) &
    direction_agree(primary_logFC, no_zone_logFC) &
    direction_agree(primary_logFC, reduced_zone_logFC) &
    direction_agree(primary_logFC, exclude_both_logFC)
)]
comparison[, stability_class := fifelse(
  sample_sensitive, "SAMPLE_SENSITIVE",
  fifelse(
    zone_dependent, "ZONE_DEPENDENT",
    fifelse(
      model_dependent, "MODEL_DEPENDENT",
      fifelse(
        primary_FDR < 0.05 & all_core_direction_stable & direction_retention_fraction >= 0.90,
        "ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY",
        fifelse(
          all_core_direction_stable & direction_retention_fraction >= 0.80,
          "DIRECTIONALLY_STABLE_LOW_PRECISION",
          "UNSTABLE"
        )
      )
    )
  )
)]
write_csv(comparison, file.path(result_root, "HD_T2_GENE_STABILITY_CLASSIFICATION.csv"))
log_msg("Gene stability classification completed.")

loo_module_dt <- rbindlist(loo_module, fill = TRUE)
module_stability <- copy(module_diff)
for (nm in names(module_camera)[-1]) {
  cam <- module_camera[[nm]]
  if (is.null(cam)) next
  tmp <- cam[, .(module_id, camera_direction, camera_p_value, camera_BH_FDR)]
  setnames(
    tmp,
    c("camera_direction", "camera_p_value", "camera_BH_FDR"),
    paste0(nm, c("_direction", "_p_value", "_BH_FDR"))
  )
  module_stability <- merge(module_stability, tmp, by = "module_id", all.x = TRUE, sort = FALSE)
}
loo_module_summary <- loo_module_dt[, .(
  loo_direction_retention = mean(camera_direction == camera_direction[which.min(camera_p_value)], na.rm = TRUE),
  loo_median_p_value = median(camera_p_value, na.rm = TRUE),
  loo_FDR_lt_0_05_fraction = mean(camera_BH_FDR < 0.05, na.rm = TRUE)
), by = module_id]
primary_dir_map <- module_diff[, .(module_id, primary_direction = camera_direction)]
loo_module_summary <- merge(loo_module_summary, primary_dir_map, by = "module_id", all.x = TRUE)
loo_module_summary <- loo_module_dt[loo_module_summary, on = "module_id", allow.cartesian = TRUE][, .(
  loo_direction_retention = mean(camera_direction == first(primary_direction), na.rm = TRUE),
  loo_median_p_value = median(camera_p_value, na.rm = TRUE),
  loo_FDR_lt_0_05_fraction = mean(camera_BH_FDR < 0.05, na.rm = TRUE)
), by = module_id]
module_stability <- merge(module_stability, loo_module_summary, by = "module_id", all.x = TRUE, sort = FALSE)
module_stability[, stable_predefined_module := (
  camera_BH_FDR < 0.05 &
    UNADJUSTED_direction == camera_direction &
    NO_ZONE_ADJUSTMENT_direction == camera_direction &
    REDUCED_ZONE_direction == camera_direction &
    EXCLUDE_BOTH_WARNING_FILES_direction == camera_direction &
    loo_direction_retention >= 0.80
)]
module_stability[, stability_label := fifelse(
  stable_predefined_module, "STABLE_FDR_SUPPORTED",
  fifelse(
    loo_direction_retention < 0.80, "SAMPLE_SENSITIVE",
    fifelse(
      NO_ZONE_ADJUSTMENT_direction != camera_direction |
        REDUCED_ZONE_direction != camera_direction,
      "ZONE_OR_MODEL_DEPENDENT",
      "NOT_FDR_SUPPORTED_OR_LOW_PRECISION"
    )
  )
)]
write_csv(module_stability, file.path(result_root, "HD_T2_MODULE_STABILITY_RESULTS.csv"))
write_csv(loo_module_dt, file.path(result_root, "06_modules", "HD_T2_MODULE_LEAVE_ONE_OUT_RESULTS.csv"))

tx_model <- fit_voom_model(
  tx_counts, meta,
  ~ disease_group + sex + BMI_centered + tissue_zone,
  "TRANSCRIPT_LEVEL_PRIMARY_SENSITIVITY"
)
if (tx_model$status != "ESTIMABLE") stop("Transcript sensitivity model failed.")
tx_res <- copy(tx_model$result)
setnames(tx_res, "feature_id", "original_feature_id")
tx_res <- merge(
  tx_res,
  tx_annotation_tested[, .(
    original_feature_id, source_enst_versioned, source_enst_stable,
    ensembl_gene_id, gene_symbol, mapped_gene_n, mapping_status
  )],
  by = "original_feature_id", all.x = TRUE, sort = FALSE
)
write_gz_csv(tx_res, file.path(result_root, "HD_T2_TRANSCRIPT_LEVEL_PRIMARY_MODEL.csv.gz"))

tx_gene <- tx_res[mapped_gene_n == 1L & !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id)]
tx_gene[, direction_nontrivial := fifelse(
  abs(log2_fold_change) < 0.10, 0L, fifelse(log2_fold_change > 0, 1L, -1L)
)]
concordance <- tx_gene[, {
  primary_gene_effect <- primary_res$log2_fold_change[match(first(ensembl_gene_id), primary_res$ensembl_gene_id)]
  nonzero_dirs <- unique(direction_nontrivial[direction_nontrivial != 0L])
  lead <- which.max(abs(moderated_t))
  .(
    tested_transcript_n = .N,
    gene_primary_log2FC = primary_gene_effect,
    transcript_same_direction_fraction = mean(
      sign(log2_fold_change) == sign(primary_gene_effect), na.rm = TRUE
    ),
    transcript_log2FC_min = min(log2_fold_change, na.rm = TRUE),
    transcript_log2FC_max = max(log2_fold_change, na.rm = TRUE),
    max_absolute_transcript_log2FC = max(abs(log2_fold_change), na.rm = TRUE),
    lead_transcript_id = original_feature_id[lead],
    lead_transcript_log2FC = log2_fold_change[lead],
    any_transcript_FDR_lt_0_05 = any(BH_FDR < 0.05, na.rm = TRUE),
    transcript_heterogeneity_warning = length(nonzero_dirs) > 1L
  )
}, by = .(ensembl_gene_id, gene_symbol)]
concordance[, concordance_label := fifelse(
  transcript_heterogeneity_warning,
  "TRANSCRIPT_HETEROGENEITY_WARNING",
  fifelse(
    transcript_same_direction_fraction >= 0.80,
    "GENE_TRANSCRIPT_DIRECTION_CONCORDANT",
    "LOW_TRANSCRIPT_CONCORDANCE"
  )
)]
write_csv(concordance, file.path(result_root, "HD_T2_GENE_TRANSCRIPT_CONCORDANCE.csv"))
write_md(
  file.path(result_root, "03_transcript_level", "HD_T2_TRANSCRIPT_SENSITIVITY_SUMMARY.md"),
  c(
    "# HD-T2 transcript sensitivity summary",
    "",
    sprintf("- ENST before filtering: %d.", nrow(raw_mat)),
    sprintf("- ENST tested: %d.", nrow(tx_counts)),
    sprintf("- Genes with tested uniquely mapped transcripts: %d.", nrow(concordance)),
    sprintf("- Genes marked TRANSCRIPT_HETEROGENEITY_WARNING: %d.", sum(concordance$transcript_heterogeneity_warning)),
    "- Transcript-level results are sensitivity evidence only and do not replace gene-level primary results.",
    "- Opposite transcript directions are not interpreted as differential splicing."
  )
)
log_msg("Transcript-level sensitivity completed.")

save_plot_pair <- function(plot, stem, width_mm = 150, height_mm = 105) {
  pdf_path <- file.path(result_root, "08_figures", paste0(stem, ".pdf"))
  png_path <- file.path(result_root, "08_figures", paste0(stem, ".png"))
  ggsave(
    pdf_path, plot = plot, device = cairo_pdf,
    width = width_mm / 25.4, height = height_mm / 25.4, units = "in"
  )
  ggsave(
    png_path, plot = plot, device = "png",
    width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
    dpi = 600, bg = "white"
  )
  data.table(
    figure_id = stem, pdf_path = pdf_path, png_path = png_path,
    width_mm = width_mm, height_mm = height_mm, dpi_png = 600
  )
}

theme_hd <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 9)
  )
palette_group <- c(
  "anal fissure surgical control" = "#6B778D",
  "hemorrhoidal disease" = "#2A7F9E"
)

figure_index <- list()
if (!is.null(model_primary$v$voom.xy) && !is.null(model_primary$v$voom.line)) {
  voom_df <- data.table(
    mean_log_count = model_primary$v$voom.xy$x,
    sqrt_residual_sd = model_primary$v$voom.xy$y
  )
  line_df <- data.table(
    mean_log_count = model_primary$v$voom.line$x,
    fitted_sqrt_residual_sd = model_primary$v$voom.line$y
  )
  p1 <- ggplot(voom_df, aes(mean_log_count, sqrt_residual_sd)) +
    geom_point(size = 0.35, alpha = 0.25, colour = "#6B778D") +
    geom_line(data = line_df, aes(y = fitted_sqrt_residual_sd), colour = "#2A7F9E", linewidth = 0.7) +
    labs(
      title = "Voom mean-variance trend",
      x = "Mean log2 expression", y = "Square-root residual SD"
    ) + theme_hd
  figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p1, "HD_T2_FIG01_VOOM_MEAN_VARIANCE")
}

qw_plot_dt <- copy(quality_weights)[order(quality_weight)]
qw_plot_dt[, patient_id := factor(patient_id, levels = patient_id)]
p2 <- ggplot(qw_plot_dt, aes(patient_id, quality_weight, fill = disease_group)) +
  geom_col(width = 0.78) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
  scale_fill_manual(values = palette_group) +
  labs(
    title = "Patient-level voom quality weights",
    x = "Patient", y = "Quality weight", fill = NULL
  ) +
  theme_hd +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p2, "HD_T2_FIG02_SAMPLE_QUALITY_WEIGHTS", 183, 105)

p3 <- ggplot(primary_res, aes(log2_fold_change)) +
  geom_histogram(bins = 80, fill = "#2A7F9E", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  labs(
    title = "Adjusted gene-effect distribution",
    x = "Adjusted log2 fold change", y = "Tested genes"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p3, "HD_T2_FIG03_ADJUSTED_EFFECT_DISTRIBUTION")

p4 <- ggplot(comparison, aes(unadjusted_logFC, primary_logFC)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#6B778D", linewidth = 0.4) +
  geom_point(alpha = 0.25, size = 0.45, colour = "#2A7F9E") +
  coord_equal() +
  labs(
    title = "Adjusted versus unadjusted effects",
    subtitle = sprintf("Spearman rho = %.3f", cor(comparison$primary_logFC, comparison$unadjusted_logFC, method = "spearman")),
    x = "Unadjusted log2 fold change", y = "Adjusted log2 fold change"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p4, "HD_T2_FIG04_ADJUSTED_VS_UNADJUSTED")

p5 <- ggplot(comparison, aes(reduced_zone_logFC, primary_logFC, colour = stability_class)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#6B778D", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.5) +
  scale_colour_manual(values = c(
    ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY = "#1B9E77",
    DIRECTIONALLY_STABLE_LOW_PRECISION = "#2A7F9E",
    ZONE_DEPENDENT = "#E28E2C",
    SAMPLE_SENSITIVE = "#D24B40",
    MODEL_DEPENDENT = "#8C6BB1",
    UNSTABLE = "#9E9E9E"
  )) +
  coord_equal() +
  labs(
    title = "Primary versus reduced-zone effects",
    x = "Reduced-zone log2 fold change",
    y = "Primary adjusted log2 fold change",
    colour = "Stability"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p5, "HD_T2_FIG05_PRIMARY_VS_REDUCED_ZONE", 183, 120)

module_plot <- copy(module_diff)
module_plot[, module_name := factor(module_name, levels = rev(module_name[order(score_log2_difference)]))]
p6 <- ggplot(module_plot, aes(score_log2_difference, module_name, size = -log10(camera_BH_FDR), colour = camera_direction)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.35) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = c(Up = "#D24B40", Down = "#3182BD")) +
  labs(
    title = "Predefined vascular-stromal module effects",
    x = "Adjusted module-score difference", y = NULL,
    size = "-log10 camera FDR", colour = "camera"
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p6, "HD_T2_FIG06_PREDEFINED_MODULE_EFFECTS", 165, 125)

stability_counts <- comparison[, .N, by = stability_class][order(-N)]
p7 <- ggplot(stability_counts, aes(reorder(stability_class, N), N, fill = stability_class)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c(
    ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY = "#1B9E77",
    DIRECTIONALLY_STABLE_LOW_PRECISION = "#2A7F9E",
    ZONE_DEPENDENT = "#E28E2C",
    SAMPLE_SENSITIVE = "#D24B40",
    MODEL_DEPENDENT = "#8C6BB1",
    UNSTABLE = "#9E9E9E"
  )) +
  labs(title = "Gene stability classification", x = NULL, y = "Tested genes") +
  theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p7, "HD_T2_FIG07_STABILITY_CLASSIFICATION")

candidate_genes <- comparison[
  stability_class %chin% c("ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY", "DIRECTIONALLY_STABLE_LOW_PRECISION")
][order(primary_FDR, -abs(primary_logFC))]
if (nrow(candidate_genes) < 6L) candidate_genes <- comparison[order(primary_FDR, -abs(primary_logFC))]
figure_robust_gene_n <- sum(
  comparison$stability_class == "ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY"
)
top_ids <- head(candidate_genes$ensembl_gene_id, 6L)
top_symbols <- gene_meta_tested$gene_symbol[match(top_ids, gene_meta_tested$ensembl_gene_id)]
top_labels <- ifelse(is.na(top_symbols) | !nzchar(top_symbols), top_ids, top_symbols)
expr_long <- data.table(
  ensembl_gene_id = rep(top_ids, times = ncol(model_primary$v$E)),
  patient_id = rep(colnames(model_primary$v$E), each = length(top_ids)),
  voom_logCPM = as.vector(model_primary$v$E[top_ids, , drop = FALSE])
)
expr_long[, gene_label := top_labels[match(ensembl_gene_id, top_ids)]]
expr_long <- merge(expr_long, meta, by = "patient_id", all.x = TRUE, sort = FALSE)
p8 <- ggplot(expr_long, aes(disease_group, voom_logCPM, colour = disease_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.35) +
  geom_jitter(width = 0.13, height = 0, size = 0.9, alpha = 0.75) +
  facet_wrap(~ gene_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = palette_group) +
  scale_x_discrete(labels = c(
    "anal fissure surgical control" = "anal fissure\nsurgical control",
    "hemorrhoidal disease" = "hemorrhoidal\ndisease"
  )) +
  labs(
    title = if (figure_robust_gene_n > 0L)
      "Patient-level expression of leading robust genes"
    else
      "Leading directionally stable, low-precision genes",
    x = NULL, y = "Voom log2 CPM", colour = NULL
  ) + theme_hd +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p8, "HD_T2_FIG08_TOP_STABLE_GENE_EXPRESSION", 183, 130)

p9 <- ggplot(expr_long, aes(disease_group, voom_logCPM, colour = disease_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.3) +
  geom_jitter(width = 0.13, height = 0, size = 0.75, alpha = 0.7) +
  facet_grid(gene_label ~ tissue_zone, scales = "free_y") +
  scale_colour_manual(values = palette_group) +
  scale_x_discrete(labels = c(
    "anal fissure surgical control" = "anal fissure\nsurgical control",
    "hemorrhoidal disease" = "hemorrhoidal\ndisease"
  )) +
  labs(
    title = if (figure_robust_gene_n > 0L)
      "Zone-stratified expression of leading robust genes"
    else
      "Zone-stratified directionally stable, low-precision genes",
    x = NULL, y = "Voom log2 CPM", colour = NULL
  ) + theme_hd +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p9, "HD_T2_FIG09_ZONE_STRATIFIED_EXPRESSION", 183, 210)

p10 <- ggplot(concordance, aes(gene_primary_log2FC, lead_transcript_log2FC, colour = concordance_label)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#6B778D", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.5) +
  scale_colour_manual(values = c(
    GENE_TRANSCRIPT_DIRECTION_CONCORDANT = "#2A7F9E",
    LOW_TRANSCRIPT_CONCORDANCE = "#E28E2C",
    TRANSCRIPT_HETEROGENEITY_WARNING = "#D24B40"
  )) +
  labs(
    title = "Gene-transcript effect concordance",
    x = "Gene-level adjusted log2 fold change",
    y = "Lead-transcript adjusted log2 fold change",
    colour = NULL
  ) + theme_hd
figure_index[[length(figure_index) + 1L]] <- save_plot_pair(p10, "HD_T2_FIG10_GENE_TRANSCRIPT_CONCORDANCE", 165, 120)

figure_index_dt <- rbindlist(figure_index, fill = TRUE)
figure_index_dt[, `:=`(
  statistical_unit = "patient/library",
  status = "HD-T2 result candidate; not frozen final manuscript figure",
  visual_QA = "PENDING_MANUAL_INSPECTION"
)]
write_csv(figure_index_dt, file.path(result_root, "08_figures", "HD_T2_FIGURE_INDEX.csv"))
write_md(
  file.path(result_root, "HD_T2_FIGURE_INDEX.md"),
  c(
    "# HD-T2 figure index",
    "",
    "- Core conclusion: patient-level mRNA effects must be shown together with tissue-zone, sample and transcript stability.",
    "- Archetype: quantitative robustness series.",
    "- Backend: R only.",
    "- Export: PDF plus 600-dpi PNG.",
    "- These are HD-T2 candidate result figures, not frozen final manuscript figures.",
    "",
    paste0("- ", figure_index_dt$figure_id, ": ", basename(figure_index_dt$pdf_path), " + ", basename(figure_index_dt$png_path))
  )
)
log_msg("Candidate figures exported in PDF and 600-dpi PNG.")

primary_sig_n <- sum(model_primary$result$BH_FDR < 0.05, na.rm = TRUE)
unadjusted_sig_n <- sum(model_unadjusted$result$BH_FDR < 0.05, na.rm = TRUE)
reduced_sig_n <- if (model_reduced$status == "ESTIMABLE") sum(model_reduced$result$BH_FDR < 0.05, na.rm = TRUE) else NA_integer_
stable_gene_n <- sum(comparison$stability_class == "ROBUST_ACROSS_PRIMARY_AND_SENSITIVITY")
directionally_stable_n <- sum(comparison$stability_class == "DIRECTIONALLY_STABLE_LOW_PRECISION")
stable_modules <- module_stability[stable_predefined_module == TRUE, module_name]
zone_dependent_fraction <- mean(comparison$stability_class == "ZONE_DEPENDENT")
sample_sensitive_fraction <- mean(comparison$stability_class == "SAMPLE_SENSITIVE")

module_supported <- length(stable_modules) > 0L
central_remodeling_modules <- c(
  "endothelial activation", "angiogenesis", "venous/vascular remodeling",
  "pericyte", "vascular smooth muscle", "fibroblast activation",
  "extracellular matrix organization", "collagen organization",
  "focal adhesion", "mechanotransduction"
)
central_module_review <- module_stability[module_name %chin% central_remodeling_modules]
central_zone_or_model_dependent_n <- sum(
  central_module_review$stability_label == "ZONE_OR_MODEL_DEPENDENT",
  na.rm = TRUE
)
central_zone_or_model_dependent_fraction <- central_zone_or_model_dependent_n /
  nrow(central_module_review)
major_zone_problem <- zone_dependent_fraction > 0.50 ||
  cor(comparison$primary_logFC, comparison$no_zone_logFC, method = "spearman", use = "complete.obs") < 0.50 ||
  central_zone_or_model_dependent_fraction >= 0.50
major_sample_problem <- sample_sensitive_fraction > 0.50

final_status <- if (major_sample_problem) {
  "HD_T2_FAIL_UNRESOLVED_CONFOUNDING"
} else if (major_zone_problem) {
  "HD_T2_PASS_ZONE_DEPENDENT_SIGNAL_REQUIRES_WORK_REVIEW"
} else if (module_supported || stable_gene_n > 0L) {
  "HD_T2_PASS_STABLE_MRNA_SIGNAL_READY_FOR_T3"
} else {
  "HD_T2_PASS_WEAK_OR_NULL_MRNA_SIGNAL_READY_FOR_T3"
}

hypothesis_support <- if (major_zone_problem) {
  paste0(
    "PARTIAL_SUPPORT_WITH_MAJOR_TISSUE_ZONE_DEPENDENCE; stable competitive modules: ",
    if (length(stable_modules)) paste(stable_modules, collapse = "; ") else "none",
    "; central remodeling modules classified zone/model-dependent: ",
    central_zone_or_model_dependent_n, "/", nrow(central_module_review)
  )
} else if (module_supported) {
  paste0("SUPPORTED_BY_STABLE_COMPETITIVE_MODULES: ", paste(stable_modules, collapse = "; "))
} else {
  "NOT_SUPPORTED_BY_FDR_STABLE_PREDEFINED_MODULES_IN_HD_T2"
}

results_ledger <- data.table(
  item = c(
    "patients_primary", "genes_aggregated", "genes_tested", "transcripts_tested",
    "primary_FDR_significant_genes", "unadjusted_FDR_significant_genes",
    "reduced_zone_FDR_significant_genes", "robust_gene_n",
    "directionally_stable_low_precision_gene_n", "stable_predefined_module_n",
    "zone_dependent_gene_fraction", "sample_sensitive_gene_fraction",
    "central_zone_or_model_dependent_module_fraction",
    "quality_weight_min", "quality_weight_max", "final_status"
  ),
  value = c(
    nrow(meta), nrow(gene_raw), nrow(gene_counts), nrow(tx_counts),
    primary_sig_n, unadjusted_sig_n, reduced_sig_n, stable_gene_n,
    directionally_stable_n, length(stable_modules),
    zone_dependent_fraction, sample_sensitive_fraction,
    central_zone_or_model_dependent_fraction,
    min(quality_weights$quality_weight), max(quality_weights$quality_weight),
    final_status
  )
)
write_csv(results_ledger, file.path(result_root, "HD_T2_RESULTS_LEDGER.md.csv"))
write_md(
  file.path(result_root, "HD_T2_RESULTS_LEDGER.md"),
  c(
    "# HD-T2 results ledger",
    "",
    paste0("- ", results_ledger$item, ": ", results_ledger$value)
  )
)

claim_boundary <- c(
  "# HD-T2 claim boundary",
  "",
  "## Allowed",
  "",
  "- Patient-level mRNA differences in hemorrhoidal disease relative to anal fissure surgical controls.",
  "- Adjusted associations after accounting for sex, continuous BMI and tissue zone.",
  "- Stability or instability across frozen models, zones and patient-level sensitivities.",
  "- Patient-level association of predefined vascular-stromal modules.",
  "",
  "## Prohibited",
  "",
  "- Healthy-versus-disease or hemorrhoid-specific claims.",
  "- Causal, diagnostic, mechanistic or cell-proportion claims.",
  "- External replication claims.",
  "- Treating transcript heterogeneity as differential splicing.",
  "- Treating low-precision within-zone significance as a primary conclusion.",
  "",
  paste0("Central vascular-stromal hypothesis in HD-T2: `", hypothesis_support, "`."),
  "HD-T3 was not started."
)
write_md(file.path(result_root, "HD_T2_CLAIM_BOUNDARY.md"), claim_boundary)

input_report <- c(
  "# HD-T2 input integrity report",
  "",
  "- Patients: 38; 20 hemorrhoidal disease and 18 anal fissure surgical controls.",
  "- Raw ENST matrix: 32,153 rows x 38 patient columns.",
  "- Column order exactly matches the formal manifest.",
  "- No NA, Inf, negative values, duplicate features or duplicate patients.",
  "- Fractional raw precision was preserved.",
  "- GSM4676457 absent tail features were not filled with zero.",
  "- All formal input SHA-256 values are stored in `01_inputs/HD_T2_INPUT_SHA256.csv`.",
  "",
  paste0("- ", integrity_dt$check, ": ", ifelse(integrity_dt$pass, "PASS", "FAIL"))
)
write_md(file.path(result_root, "HD_T2_INPUT_INTEGRITY_REPORT.md"), input_report)

run_end <- Sys.time()
pkg_names <- c("limma", "edgeR", "statmod", "data.table", "ggplot2", "digest")
pkg_versions <- data.table(
  package = pkg_names,
  version = vapply(pkg_names, function(p) as.character(packageVersion(p)), character(1))
)
write_csv(pkg_versions, file.path(result_root, "09_environment", "HD_T2_PACKAGE_VERSIONS.csv"))
capture.output(sessionInfo(), file = file.path(result_root, "09_environment", "HD_T2_SESSION_INFO.txt"))
env_lock <- c(
  "# HD-T2 environment lock",
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
  "- HD-T1 matrices and source files remained read-only."
)
write_md(file.path(result_root, "HD_T2_ENVIRONMENT_LOCK.md"), env_lock)

executive_summary <- c(
  "# HD-T2 executive summary",
  "",
  paste0("**Final state:** `", final_status, "`"),
  "",
  sprintf("- Primary analysis patients: %d.", nrow(meta)),
  sprintf("- Aggregated genes: %d; tested genes: %d.", nrow(gene_raw), nrow(gene_counts)),
  "- Primary model: `~ disease_group + sex + BMI_centered + tissue_zone`.",
  sprintf("- Quality-weight range: %.4f to %.4f.", min(quality_weights$quality_weight), max(quality_weights$quality_weight)),
  sprintf("- Adjusted-model FDR-significant genes: %d.", primary_sig_n),
  sprintf("- Unadjusted-model FDR-significant genes: %d.", unadjusted_sig_n),
  sprintf("- Reduced-zone FDR-significant genes: %s.", as.character(reduced_sig_n)),
  sprintf("- Robust genes: %d.", stable_gene_n),
  sprintf("- Directionally stable low-precision genes: %d.", directionally_stable_n),
  sprintf("- Stable predefined competitive modules: %s.", if (length(stable_modules)) paste(stable_modules, collapse = "; ") else "none"),
  sprintf("- Zone-dependent gene fraction: %.3f.", zone_dependent_fraction),
  sprintf("- Sample-sensitive gene fraction: %.3f.", sample_sensitive_fraction),
  sprintf("- Central remodeling modules classified zone/model-dependent: %d/%d.", central_zone_or_model_dependent_n, nrow(central_module_review)),
  paste0("- Central hypothesis: ", hypothesis_support, "."),
  "- Formal miRNA differential analysis: NOT RUN.",
  "- HD-T3 automatic start: NO."
)
write_md(file.path(result_root, "HD_T2_EXECUTIVE_SUMMARY.md"), executive_summary)

final_gate <- c(
  "# HD-T2 final gate",
  "",
  "## Final status",
  "",
  paste0("`", final_status, "`"),
  "",
  sprintf("- Final primary-analysis patients: %d.", nrow(meta)),
  sprintf("- Final tested genes: %d.", nrow(gene_counts)),
  "- Primary adjusted model: `~ disease_group + sex + BMI_centered + tissue_zone`.",
  sprintf("- Quality weight range: %.6f to %.6f.", min(quality_weights$quality_weight), max(quality_weights$quality_weight)),
  sprintf("- Adjusted-model FDR-significant genes: %d.", primary_sig_n),
  sprintf("- Unadjusted-model FDR-significant genes: %d.", unadjusted_sig_n),
  sprintf("- Reduced-zone FDR-significant genes: %s.", as.character(reduced_sig_n)),
  sprintf("- Robust stable genes: %d.", stable_gene_n),
  sprintf("- Stable predefined competitive modules: %s.", if (length(stable_modules)) paste(stable_modules, collapse = "; ") else "none"),
  sprintf("- Tissue-zone impact: zone-dependent gene fraction %.3f.", zone_dependent_fraction),
  sprintf("- Central remodeling modules classified zone/model-dependent: %d/%d.", central_zone_or_model_dependent_n, nrow(central_module_review)),
  sprintf("- Warning-sample impact: sample-sensitive gene fraction %.3f.", sample_sensitive_fraction),
  paste0("- Vascular-stromal remodeling hypothesis: ", hypothesis_support, "."),
  "- HD-T3 may be considered only after user review of this gate.",
  "- HD-T3 automatically started: NO."
)
write_md(file.path(result_root, "HD_T2_FINAL_GATE.md"), final_gate)
file.copy(
  file.path(result_root, "HD_T2_FINAL_GATE.md"),
  file.path(result_root, "10_gate", "HD_T2_FINAL_GATE.md"),
  overwrite = TRUE
)

log_msg("HD-T2 completed with status ", final_status, ". HD-T3 not started.")
file.copy(log_file, file.path(result_root, "HD_T2_ANALYSIS_LOG.md"), overwrite = TRUE)

all_outputs <- list.files(result_root, recursive = TRUE, full.names = TRUE)
all_outputs <- all_outputs[file.info(all_outputs)$isdir == FALSE]
hash_index <- data.table(
  relative_path = substring(normalizePath(all_outputs, winslash = "/", mustWork = TRUE), nchar(project_root) + 2L),
  file_size_bytes = as.numeric(file.info(all_outputs)$size),
  sha256 = vapply(all_outputs, sha256, character(1))
)
hash_index <- hash_index[relative_path != "HD_T2_SHA256_INDEX.csv"]
write_csv(hash_index, file.path(result_root, "HD_T2_SHA256_INDEX.csv"))

cat(final_status, "\n")
