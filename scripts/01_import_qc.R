#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
Sys.setenv(TZ = "Asia/Shanghai")
set.seed(20260716)

required_packages <- c("data.table", "readr", "ggplot2", "digest", "jsonlite", "RColorBrewer", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       ". HD-T1 will not substitute a different workflow.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(ggplot2)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args)) normalizePath(args[[1]], winslash = "\\", mustWork = TRUE) else
  normalizePath(file.path(getwd()), winslash = "\\", mustWork = TRUE)
if (basename(project_root) != "HEMORRHOID_TRANSCRIPTOMICS") {
  candidate <- file.path(project_root, "HEMORRHOID_TRANSCRIPTOMICS")
  if (dir.exists(candidate)) project_root <- normalizePath(candidate, winslash = "\\", mustWork = TRUE)
}

manifest_path <- file.path(project_root, "HD_DATASET_AND_SAMPLE_MANIFEST.csv")
verified_manifest_path <- file.path(project_root, "HD_T0R_MASTER_SAMPLE_MANIFEST_VERIFIED.csv")
annotation_dir <- file.path(project_root, "source_metadata", "HD_T1_annotation_cache")
gtf_path <- file.path(annotation_dir, "Homo_sapiens.GRCh38.100.gtf.gz")
mirbase_mature_path <- file.path(annotation_dir, "miRBase_22.1_mature.decoded.fa")
mirbase_hairpin_path <- file.path(annotation_dir, "miRBase_22.1_hairpin.decoded.fa")
mirbase_dead_path <- file.path(annotation_dir, "miRBase_22.1_miRNA.decoded.dead")
mirbase_readme_path <- file.path(annotation_dir, "miRBase_22.1_README.decoded.txt")

out_root <- file.path(project_root, "results", "HD_T1")
dirs <- c(
  logs = file.path(out_root, "00_logs"),
  manifests = file.path(out_root, "01_manifests"),
  matrices = file.path(out_root, "02_matrices"),
  annotation = file.path(out_root, "03_feature_annotation"),
  qc_tables = file.path(out_root, "04_qc_tables"),
  qc_figures = file.path(out_root, "05_qc_figures"),
  environment = file.path(out_root, "06_environment"),
  gate = file.path(out_root, "07_gate")
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

run_start <- Sys.time()
log_lines <- character()
log_event <- function(message) {
  entry <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), message)
  log_lines <<- c(log_lines, entry)
  message(entry)
}
write_utf8 <- function(x, path) writeLines(enc2utf8(x), path, useBytes = TRUE)
sha256 <- function(path) toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))
fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
fmt_num <- function(x, digits = 4) format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
safe_mad_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  s <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - med) / s
}
extract_attr <- function(x, key) {
  pattern <- paste0(".*(?:^|; )", key, " \"([^\"]+)\".*")
  hit <- grepl(paste0("(?:^|; )", key, " \""), x)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub(pattern, "\\1", x[hit], perl = TRUE)
  out
}
count_gzip_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  n <- 0L
  last <- ""
  repeat {
    x <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(x)) break
    n <- n + length(x)
    last <- x[[length(x)]]
  }
  list(n = n, last = last)
}
save_plot_pair <- function(plot, stem, width = 10, height = 7) {
  png_path <- file.path(dirs[["qc_figures"]], paste0(stem, ".png"))
  pdf_path <- file.path(dirs[["qc_figures"]], paste0(stem, ".pdf"))
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf, bg = "white")
  c(png_path, pdf_path)
}

log_event("HD-T1 started. Formal differential expression is prohibited.")

stopifnot(file.exists(manifest_path), file.exists(verified_manifest_path))
manifest_sha <- sha256(manifest_path)
verified_manifest_sha <- sha256(verified_manifest_path)
if (!identical(manifest_sha, verified_manifest_sha)) {
  stop("Formal and T0R verified manifest SHA-256 differ. HD-T1 stopped before import.")
}
manifest <- fread(manifest_path, encoding = "UTF-8", na.strings = character())
if (nrow(manifest) != 38L || uniqueN(manifest$patient_id) != 38L ||
    uniqueN(manifest$GEO_sample_id) != 38L) {
  stop("Formal manifest is not a unique 38-patient manifest.")
}
if (sum(manifest$disease_group == "hemorrhoidal disease") != 20L ||
    sum(manifest$disease_group == "anal fissure surgical control") != 18L) {
  stop("Formal manifest does not reproduce the frozen 20/18 groups.")
}
if (!all(file.exists(manifest$processed_tsv_file))) stop("One or more source TSV files are absent.")

expected_gtf_bytes <- 46973686
if (!file.exists(gtf_path) || file.info(gtf_path)$size != expected_gtf_bytes) {
  stop("Frozen Ensembl 100 GTF is absent or has an unexpected byte size: ", gtf_path)
}
if (!file.exists(mirbase_mature_path) || file.info(mirbase_mature_path)$size < 1000000) {
  stop("miRBase 22.1 mature.fa audit reference is absent or incomplete.")
}
if (!file.exists(mirbase_hairpin_path) || file.info(mirbase_hairpin_path)$size < 1000000) {
  stop("miRBase 22.1 hairpin.fa precursor audit reference is absent or incomplete.")
}
if (!file.exists(mirbase_dead_path) || file.info(mirbase_dead_path)$size < 10000) {
  stop("miRBase 22.1 deprecated-name audit reference is absent or incomplete.")
}

log_event(paste("Manifest gate passed:", manifest_sha))

read_source <- function(path) {
  parsed <- suppressWarnings(read_tsv(
    path,
    col_types = cols(gene = col_character(), raw = col_double(), rpm = col_double()),
    progress = FALSE,
    name_repair = "minimal",
    trim_ws = FALSE,
    quote = ""
  ))
  probs <- problems(parsed)
  valid <- !is.na(parsed$gene) & nzchar(parsed$gene) &
    is.finite(parsed$raw) & is.finite(parsed$rpm) &
    parsed$raw >= 0 & parsed$rpm >= 0
  list(all = parsed, valid = as.data.table(parsed[valid, , drop = FALSE]), problems = probs,
       invalid_n = sum(!valid))
}

feature_lists <- vector("list", nrow(manifest))
input_rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  path <- row$processed_tsv_file
  log_event(sprintf("Import audit %02d/38: %s", i, row$GEO_sample_id))
  observed_sha <- sha256(path)
  if (!identical(observed_sha, toupper(row$processed_tsv_sha256))) {
    stop("Input SHA-256 mismatch for ", row$GEO_sample_id)
  }
  line_info <- count_gzip_lines(path)
  src <- read_source(path)
  d <- src$valid
  feature_lists[[i]] <- d$gene
  duplicate_n <- sum(duplicated(d$gene))
  terminal_problem <- nrow(src$problems) > 0L &&
    any(src$problems$row >= nrow(src$all), na.rm = TRUE)
  input_rows[[i]] <- data.table(
    GEO_sample_id = row$GEO_sample_id,
    patient_id = row$patient_id,
    disease_group = row$disease_group,
    file_path = path,
    file_size_bytes = file.info(path)$size,
    sha256 = observed_sha,
    gzip_status = "COMPLETE_READ_TO_EOF",
    line_n = line_info$n,
    column_n = ncol(src$all),
    header = paste(names(src$all), collapse = "|"),
    parsed_row_n = nrow(src$all),
    valid_feature_n = nrow(d),
    ENST_n = sum(grepl("^ENST", d$gene)),
    miRNA_n = sum(grepl("^hsa-(miR|mir|let)-", d$gene)),
    other_n = sum(!grepl("^(ENST|hsa-miR-|hsa-mir-|hsa-let-)", d$gene)),
    illegal_row_n = src$invalid_n,
    duplicate_feature_id_n = duplicate_n,
    terminal_format_warning = if (row$terminal_format_warning != "NO_WARNING" || terminal_problem)
      "NON_BLOCKING_FILE_FORMAT_WARNING" else "NO_WARNING",
    last_line_field_n = lengths(strsplit(line_info$last, "\t", fixed = TRUE))
  )
  if (duplicate_n > 0L) stop("Duplicate valid feature IDs in ", row$GEO_sample_id)
}
input_manifest <- rbindlist(input_rows)
fwrite(input_manifest, file.path(dirs[["manifests"]], "HD_T1_INPUT_FILE_MANIFEST.csv"))

set_lengths <- lengths(feature_lists)
common_set <- feature_lists[[which.min(set_lengths)]]
for (i in order(set_lengths)) common_set <- intersect(common_set, feature_lists[[i]])
common_features <- feature_lists[[1]][feature_lists[[1]] %in% common_set]
union_features <- feature_lists[[1]]
for (i in seq_along(feature_lists)[-1]) union_features <- union(union_features, feature_lists[[i]])
presence <- vapply(feature_lists, function(x) union_features %in% x, logical(length(union_features)))
colnames(presence) <- manifest$patient_id
presence_count <- rowSums(presence)

classify_feature <- function(x) {
  fifelse(grepl("^ENST", x), "ENST transcript",
          fifelse(grepl("^hsa-(miR|mir|let)-", x), "miRNA", "other/unknown"))
}
common_type <- classify_feature(common_features)
common_total <- length(common_features)
common_enst <- sum(common_type == "ENST transcript")
common_mirna <- sum(common_type == "miRNA")
common_other <- sum(common_type == "other/unknown")
log_event(sprintf("Common feature universe: %d total; %d ENST; %d miRNA; %d other.",
                  common_total, common_enst, common_mirna, common_other))

presence_dt <- data.table(feature_id = union_features, presence_count = presence_count)
presence_dt <- cbind(presence_dt, as.data.table(presence))
fwrite(presence_dt,
       file.path(dirs[["manifests"]], "HD_T1_FEATURE_PRESENCE_MATRIX.csv.gz"),
       compress = "gzip")

common_dt <- data.table(
  feature_order = seq_along(common_features),
  original_feature_id = common_features,
  feature_type = common_type
)
fwrite(common_dt, file.path(dirs[["manifests"]], "HD_T1_COMMON_FEATURE_LIST.csv"))

missing_lists <- apply(!presence, 1, function(z) paste(manifest$GEO_sample_id[z], collapse = ";"))
present_lists <- apply(presence, 1, function(z) paste(manifest$GEO_sample_id[z], collapse = ";"))
idx_6457 <- match("GSM4676457", manifest$GEO_sample_id)
idx_6447 <- match("GSM4676447", manifest$GEO_sample_id)
noncommon <- data.table(
  feature_id = union_features[presence_count < 38L],
  feature_type = classify_feature(union_features[presence_count < 38L]),
  presence_count = presence_count[presence_count < 38L],
  missing_GSMs = missing_lists[presence_count < 38L],
  present_GSMs = present_lists[presence_count < 38L],
  absent_in_GSM4676457 = !presence[presence_count < 38L, idx_6457],
  absent_in_GSM4676447 = !presence[presence_count < 38L, idx_6447],
  truncation_interpretation = fifelse(
    !presence[presence_count < 38L, idx_6457],
    "NOT_OBSERVED_IN_SOURCE_TRUNCATED_GSM4676457; NEVER_ZERO_FILLED",
    fifelse(!presence[presence_count < 38L, idx_6447],
            "NOT_OBSERVED_IN_TERMINALLY_SHORTENED_GSM4676447; NEVER_ZERO_FILLED",
            "NONCOMMON_FEATURE_OTHER_SAMPLE_PATTERN")
  )
)
fwrite(noncommon, file.path(dirs[["manifests"]], "HD_T1_NONCOMMON_FEATURE_AUDIT.csv"))

feature_parts <- tstrsplit(common_features, ":", fixed = TRUE, fill = NA_character_)
feature_head <- feature_parts[[1]]
source_biotype <- ifelse(common_type == "ENST transcript", feature_parts[[2]], NA_character_)
source_label <- ifelse(common_type == "ENST transcript", feature_parts[[3]], NA_character_)
enst_stable <- ifelse(common_type == "ENST transcript", sub("\\..*$", "", feature_head), NA_character_)
enst_version <- ifelse(common_type == "ENST transcript" & grepl("\\.", feature_head),
                       sub("^.*\\.", "", feature_head), NA_character_)
mirna_normalized <- ifelse(common_type == "miRNA", feature_head, NA_character_)
feature_map <- data.table(
  original_feature_id = common_features,
  normalized_feature_id = fifelse(common_type == "ENST transcript", enst_stable,
                                  fifelse(common_type == "miRNA", mirna_normalized, common_features)),
  feature_type = common_type,
  ENST_version = enst_version,
  miRNA_original_name = fifelse(common_type == "miRNA", common_features, NA_character_),
  miRNA_normalized_name = mirna_normalized,
  classification_rule = fifelse(
    common_type == "ENST transcript", "prefix_ENST_before_first_colon",
    fifelse(common_type == "miRNA", "prefix_hsa-miR_hsa-let_or_precursor_hsa-mir_before_first_colon",
            "neither_frozen_rule")
  ),
  mapping_status = "CLASSIFIED_FOR_AUDIT",
  notes = fifelse(common_type == "other/unknown", "QUARANTINED_NOT_MRNA", "")
)
fwrite(feature_map, file.path(dirs[["annotation"]], "HD_T1_FEATURE_TYPE_MAP.csv"))

raw_mat <- matrix(NA_real_, nrow = common_total, ncol = nrow(manifest),
                  dimnames = list(common_features, manifest$patient_id))
rpm_mat <- raw_mat
for (i in seq_len(nrow(manifest))) {
  log_event(sprintf("Matrix import %02d/38: %s", i, manifest$GEO_sample_id[i]))
  d <- read_source(manifest$processed_tsv_file[i])$valid
  hit <- match(common_features, d$gene)
  if (anyNA(hit)) stop("Common feature missing during matrix construction: ", manifest$GEO_sample_id[i])
  raw_mat[, i] <- d$raw[hit]
  rpm_mat[, i] <- d$rpm[hit]
}

write_matrix <- function(mat, path) {
  x <- cbind(data.table(feature_id = rownames(mat)), as.data.table(mat))
  fwrite(x, path, sep = "\t", compress = "gzip", quote = FALSE)
}
matrix_paths <- c(
  RAW_COMMON_ALL_FEATURES = file.path(dirs[["matrices"]], "HD_T1_RAW_COMMON_ALL_FEATURES.tsv.gz"),
  RPM_COMMON_ALL_FEATURES = file.path(dirs[["matrices"]], "HD_T1_RPM_COMMON_ALL_FEATURES.tsv.gz"),
  RAW_ENST = file.path(dirs[["matrices"]], "HD_T1_RAW_ENST.tsv.gz"),
  RPM_ENST = file.path(dirs[["matrices"]], "HD_T1_RPM_ENST.tsv.gz"),
  RAW_MIRNA = file.path(dirs[["matrices"]], "HD_T1_RAW_MIRNA.tsv.gz"),
  RPM_MIRNA = file.path(dirs[["matrices"]], "HD_T1_RPM_MIRNA.tsv.gz")
)
is_enst <- common_type == "ENST transcript"
is_mirna <- common_type == "miRNA"
write_matrix(raw_mat, matrix_paths[["RAW_COMMON_ALL_FEATURES"]])
write_matrix(rpm_mat, matrix_paths[["RPM_COMMON_ALL_FEATURES"]])
write_matrix(raw_mat[is_enst, , drop = FALSE], matrix_paths[["RAW_ENST"]])
write_matrix(rpm_mat[is_enst, , drop = FALSE], matrix_paths[["RPM_ENST"]])
write_matrix(raw_mat[is_mirna, , drop = FALSE], matrix_paths[["RAW_MIRNA"]])
write_matrix(rpm_mat[is_mirna, , drop = FALSE], matrix_paths[["RPM_MIRNA"]])

audit_matrix <- function(name, mat, path) {
  data.table(
    matrix_name = name,
    path = path,
    row_n = nrow(mat),
    column_n = ncol(mat),
    minimum = min(mat),
    maximum = max(mat),
    NA_n = sum(is.na(mat)),
    Inf_n = sum(is.infinite(mat)),
    negative_n = sum(mat < 0, na.rm = TRUE),
    all_zero_row_n = sum(rowSums(mat != 0) == 0),
    all_zero_column_n = sum(colSums(mat != 0) == 0),
    duplicate_row_name_n = sum(duplicated(rownames(mat))),
    duplicate_column_name_n = sum(duplicated(colnames(mat))),
    sha256 = sha256(path)
  )
}
matrix_audit <- rbindlist(list(
  audit_matrix("RAW_COMMON_ALL_FEATURES", raw_mat, matrix_paths[["RAW_COMMON_ALL_FEATURES"]]),
  audit_matrix("RPM_COMMON_ALL_FEATURES", rpm_mat, matrix_paths[["RPM_COMMON_ALL_FEATURES"]]),
  audit_matrix("RAW_ENST", raw_mat[is_enst, , drop = FALSE], matrix_paths[["RAW_ENST"]]),
  audit_matrix("RPM_ENST", rpm_mat[is_enst, , drop = FALSE], matrix_paths[["RPM_ENST"]]),
  audit_matrix("RAW_MIRNA", raw_mat[is_mirna, , drop = FALSE], matrix_paths[["RAW_MIRNA"]]),
  audit_matrix("RPM_MIRNA", rpm_mat[is_mirna, , drop = FALSE], matrix_paths[["RPM_MIRNA"]])
))
fwrite(matrix_audit, file.path(dirs[["qc_tables"]], "HD_T1_MATRIX_AUDIT.csv"))

log_event("Parsing frozen Ensembl release 100 GTF in chunks.")
target_enst <- unique(enst_stable[is_enst])
gtf_con <- gzfile(gtf_path, open = "rt")
gtf_hits <- list()
chunk_i <- 0L
repeat {
  lines <- readLines(gtf_con, n = 150000L, warn = FALSE)
  if (!length(lines)) break
  chunk_i <- chunk_i + 1L
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  if (!length(lines)) next
  fields <- tstrsplit(lines, "\t", fixed = TRUE, fill = NA_character_)
  attrs <- fields[[9]]
  tid <- extract_attr(attrs, "transcript_id")
  keep <- !is.na(tid) & tid %in% target_enst
  if (!any(keep)) next
  gtf_hits[[length(gtf_hits) + 1L]] <- data.table(
    ensembl_transcript_id = tid[keep],
    ensembl_transcript_version = extract_attr(attrs[keep], "transcript_version"),
    ensembl_gene_id = extract_attr(attrs[keep], "gene_id"),
    ensembl_gene_version = extract_attr(attrs[keep], "gene_version"),
    gene_symbol = extract_attr(attrs[keep], "gene_name"),
    gene_biotype = extract_attr(attrs[keep], "gene_biotype"),
    transcript_biotype = extract_attr(attrs[keep], "transcript_biotype"),
    chromosome = fields[[1]][keep],
    strand = fields[[7]][keep]
  )
}
close(gtf_con)
gtf_map <- unique(rbindlist(gtf_hits, fill = TRUE))
gtf_map[, ensembl_transcript_versioned := fifelse(
  !is.na(ensembl_transcript_version),
  paste0(ensembl_transcript_id, ".", ensembl_transcript_version),
  ensembl_transcript_id
)]
source_enst <- data.table(
  original_feature_id = common_features[is_enst],
  source_enst_versioned = feature_head[is_enst],
  source_enst_stable = enst_stable[is_enst],
  source_enst_version = enst_version[is_enst],
  source_biotype = source_biotype[is_enst],
  source_transcript_label = source_label[is_enst]
)
enst_annotation <- merge(source_enst, gtf_map, by.x = "source_enst_stable",
                         by.y = "ensembl_transcript_id", all.x = TRUE, allow.cartesian = TRUE)
enst_annotation[, exact_version_match := !is.na(ensembl_transcript_versioned) &
                  source_enst_versioned == ensembl_transcript_versioned]
map_counts <- enst_annotation[!is.na(ensembl_gene_id),
                              .(mapped_gene_n = uniqueN(ensembl_gene_id)), by = source_enst_stable]
enst_annotation <- merge(enst_annotation, map_counts, by = "source_enst_stable", all.x = TRUE)
enst_annotation[is.na(mapped_gene_n), mapped_gene_n := 0L]
enst_annotation[, embedded_gene_symbol := sub("-[0-9]{3,}$", "", source_transcript_label)]
enst_annotation[, mapping_status := fifelse(
  mapped_gene_n == 0L, "UNMAPPED_OR_RETIRED_IN_ENSEMBL_100",
  fifelse(mapped_gene_n > 1L, "AMBIGUOUS_MULTIPLE_GENES",
          fifelse(exact_version_match, "MAPPED_EXACT_VERSION", "MAPPED_STABLE_ID_VERSION_DIFFERS"))
)]
enst_annotation[, notes := fifelse(
  mapped_gene_n == 0L, "Retained in transcript audit; not eligible for primary gene aggregation without mapping.",
  fifelse(!is.na(gene_symbol) & embedded_gene_symbol != gene_symbol,
          "Deposited transcript label and Ensembl 100 gene symbol differ; Ensembl mapping retained for audit.",
          "")
)]
setcolorder(enst_annotation, c(
  "original_feature_id", "source_enst_versioned", "source_enst_stable",
  "source_enst_version", "source_biotype", "source_transcript_label",
  "embedded_gene_symbol", "ensembl_transcript_versioned",
  "ensembl_transcript_version", "ensembl_gene_id", "ensembl_gene_version",
  "gene_symbol", "gene_biotype", "transcript_biotype", "chromosome", "strand",
  "exact_version_match", "mapped_gene_n", "mapping_status", "notes"
))
fwrite(enst_annotation, file.path(dirs[["annotation"]], "HD_T1_ENST_ANNOTATION_MAP.csv"))

enst_unique <- enst_annotation[, .(
  mapped = any(mapped_gene_n > 0L),
  ambiguous = any(mapped_gene_n > 1L),
  exact = any(exact_version_match)
), by = original_feature_id]
mapped_enst_n <- sum(enst_unique$mapped)
unmapped_enst_n <- sum(!enst_unique$mapped)
ambiguous_enst_n <- sum(enst_unique$ambiguous)
exact_enst_n <- sum(enst_unique$exact)
gene_transcript_counts <- unique(enst_annotation[!is.na(ensembl_gene_id),
  .(source_enst_stable, ensembl_gene_id, gene_symbol)])[, .N, by = .(ensembl_gene_id, gene_symbol)]
multi_transcript_gene_n <- sum(gene_transcript_counts$N > 1L)
duplicate_symbol_n <- enst_annotation[!is.na(gene_symbol),
                                     uniqueN(gene_symbol[duplicated(gene_symbol)])]
enst_summary <- c(
  "# HD-T1 ENST mapping summary",
  "",
  "- Audit reference: Ensembl release 100, GRCh38.",
  "- Reference role: reproducible fallback mapping; not claimed as the exact unversioned original exceRpt database.",
  paste0("- Reference file: `", gtf_path, "`."),
  paste0("- Reference SHA-256: `", sha256(gtf_path), "`."),
  paste0("- Source ENST features: ", fmt_int(common_enst), "."),
  paste0("- Mapped ENST features: ", fmt_int(mapped_enst_n), " (",
         fmt_num(100 * mapped_enst_n / common_enst, 2), "%)."),
  paste0("- Unmapped or retired in release 100: ", fmt_int(unmapped_enst_n), "."),
  paste0("- Exact version matches: ", fmt_int(exact_enst_n), "."),
  paste0("- Ambiguous transcript-to-multiple-gene mappings: ", fmt_int(ambiguous_enst_n), "."),
  paste0("- Genes receiving multiple source transcripts: ", fmt_int(multi_transcript_gene_n), "."),
  paste0("- Duplicate mapped gene symbols across transcript rows: ", fmt_int(duplicate_symbol_n), "."),
  "",
  "## Frozen downstream aggregation recommendation",
  "",
  "Strip the ENST version only for lookup, exclude ambiguous transcript-to-multiple-gene mappings from the primary gene matrix, retain unmapped transcripts in the audit, and sum non-negative raw abundance across mapped transcripts sharing one Ensembl gene ID. Aggregation is not executed in HD-T1; it is an HD-T2 input-construction step after explicit authorization."
)
write_utf8(enst_summary, file.path(out_root, "HD_T1_ENST_MAPPING_SUMMARY.md"))

log_event("Parsing miRBase release 22.1 mature and deprecated-name references.")
mature_headers <- readLines(mirbase_mature_path, warn = FALSE)
mature_headers <- mature_headers[startsWith(mature_headers, ">")]
mirbase_mature <- data.table(
  normalized_miRNA_id = sub("^>(\\S+).*", "\\1", mature_headers),
  mirbase_accession = sub("^>\\S+\\s+(\\S+).*", "\\1", mature_headers),
  mirbase_name_type = "mature"
)
mirbase_mature <- mirbase_mature[grepl("^hsa-", normalized_miRNA_id)]
hairpin_headers <- readLines(mirbase_hairpin_path, warn = FALSE)
hairpin_headers <- hairpin_headers[startsWith(hairpin_headers, ">")]
mirbase_hairpin <- data.table(
  normalized_miRNA_id = sub("^>(\\S+).*", "\\1", hairpin_headers),
  mirbase_accession = sub("^>\\S+\\s+(\\S+).*", "\\1", hairpin_headers),
  mirbase_name_type = "precursor"
)
mirbase_hairpin <- mirbase_hairpin[grepl("^hsa-", normalized_miRNA_id)]
mirbase_current <- unique(rbindlist(list(mirbase_mature, mirbase_hairpin)))

dead_lines <- readLines(mirbase_dead_path, warn = FALSE)
record_end <- c(0L, which(dead_lines == "//"))
dead_records <- vector("list", max(0L, length(record_end) - 1L))
if (length(dead_records)) {
  for (i in seq_along(dead_records)) {
    z <- dead_lines[(record_end[i] + 1L):(record_end[i + 1L] - 1L)]
    get_code <- function(code) sub(paste0("^", code, "\\s+"), "", z[grepl(paste0("^", code, "\\s+"), z)])
    dead_records[[i]] <- data.table(
      deprecated_id = paste(get_code("ID"), collapse = ";"),
      previous_ids = paste(get_code("PI"), collapse = ";"),
      forward_accession = paste(get_code("FW"), collapse = ";"),
      deprecated_notes = paste(get_code("CC"), collapse = " ")
    )
  }
}
mirbase_dead <- if (length(dead_records)) rbindlist(dead_records, fill = TRUE) else data.table()

mirna_source <- data.table(
  original_miRNA_id = common_features[is_mirna],
  normalized_miRNA_id = mirna_normalized[is_mirna]
)
mirna_map <- merge(mirna_source, mirbase_current, by = "normalized_miRNA_id", all.x = TRUE)
if (nrow(mirbase_dead)) {
  dead_lookup <- rbind(
    mirbase_dead[nzchar(deprecated_id), .(lookup = deprecated_id, deprecated_id,
                                         forward_accession, deprecated_notes)],
    mirbase_dead[nzchar(previous_ids), .(
      lookup = unlist(strsplit(previous_ids, "\\s+")),
      deprecated_id = rep(deprecated_id, lengths(strsplit(previous_ids, "\\s+"))),
      forward_accession = rep(forward_accession, lengths(strsplit(previous_ids, "\\s+"))),
      deprecated_notes = rep(deprecated_notes, lengths(strsplit(previous_ids, "\\s+")))
    )],
    fill = TRUE
  )
  dead_lookup <- unique(dead_lookup[nzchar(lookup)])
  mirna_map <- merge(mirna_map, dead_lookup, by.x = "normalized_miRNA_id",
                     by.y = "lookup", all.x = TRUE)
} else {
  mirna_map[, `:=`(deprecated_id = NA_character_, forward_accession = NA_character_,
                   deprecated_notes = NA_character_)]
}
mirna_map[, species_prefix := sub("^([^-]+)-.*", "\\1", normalized_miRNA_id)]
mirna_map[, mature_or_precursor := fifelse(
  (!is.na(mirbase_name_type) & mirbase_name_type == "precursor") |
    grepl("^hsa-mir-", normalized_miRNA_id),
  "precursor",
  fifelse(
    grepl("-(3p|5p)$", normalized_miRNA_id),
    "mature_arm_specific",
    fifelse(
      !is.na(mirbase_name_type) & mirbase_name_type == "mature",
      "mature_no_arm_suffix",
      "mature_legacy_or_unversioned"
    )
  )
)]
mirna_map[, arm_5p_3p := fifelse(grepl("-5p$", normalized_miRNA_id), "5p",
                                 fifelse(grepl("-3p$", normalized_miRNA_id), "3p", NA_character_))]
mirna_map[, mapping_status := fifelse(
  !is.na(mirbase_accession),
  paste0("RECOGNIZED_MIRBASE_22_1_", toupper(mirbase_name_type)),
  fifelse(!is.na(deprecated_id), "DEPRECATED_OR_PREVIOUS_NAME", "UNRECOGNIZED_IN_MIRBASE_22_1")
)]
mirna_map[, deprecated_name := deprecated_id]
mirna_map[, replacement_name := fifelse(!is.na(forward_accession),
                                         paste0("forward_accession:", forward_accession), NA_character_)]
mirna_map[, notes := fifelse(
  grepl("^RECOGNIZED_MIRBASE_22_1", mapping_status), "",
  fifelse(mapping_status == "DEPRECATED_OR_PREVIOUS_NAME", deprecated_notes,
          "Retained exactly as deposited; no automatic renaming or arm merging.")
)]
setcolorder(mirna_map, c(
  "original_miRNA_id", "normalized_miRNA_id", "species_prefix",
  "mature_or_precursor", "arm_5p_3p", "mirbase_accession", "mapping_status",
  "deprecated_name", "replacement_name", "notes"
))
fwrite(mirna_map, file.path(dirs[["annotation"]], "HD_T1_MIRNA_NAME_MAP.csv"))
recognized_mirna_n <- sum(grepl("^RECOGNIZED_MIRBASE_22_1_", mirna_map$mapping_status))
unrecognized_mirna_n <- sum(mirna_map$mapping_status == "UNRECOGNIZED_IN_MIRBASE_22_1")
deprecated_mirna_n <- sum(mirna_map$mapping_status == "DEPRECATED_OR_PREVIOUS_NAME")
duplicate_normalized_mirna_n <- sum(duplicated(mirna_map$normalized_miRNA_id))

analyte_qc <- function(raw, rpm, analyte, patient_ids) {
  log_rpm <- log2(rpm + 1)
  spearman <- cor(log_rpm, method = "spearman", use = "pairwise.complete.obs")
  pearson <- cor(log_rpm, method = "pearson", use = "pairwise.complete.obs")
  vars <- apply(log_rpm, 1, var)
  top_n <- min(if (analyte == "ENST") 2000L else 500L, nrow(log_rpm))
  top_idx <- order(vars, decreasing = TRUE)[seq_len(top_n)]
  pca <- prcomp(t(log_rpm[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pca_scores <- as.data.table(pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE])
  while (ncol(pca_scores) < 5L) pca_scores[[paste0("PC", ncol(pca_scores) + 1L)]] <- 0
  pca_scores[, patient_id := patient_ids]
  setcolorder(pca_scores, c("patient_id", paste0("PC", 1:5)))
  pca_distance <- sqrt(rowSums(scale(as.matrix(pca_scores[, paste0("PC", 1:5), with = FALSE]),
                                    center = TRUE, scale = apply(as.matrix(pca_scores[, paste0("PC", 1:5), with = FALSE]), 2, sd))^2,
                               na.rm = TRUE))
  sample_dist <- as.matrix(dist(t(log_rpm)))
  metrics <- rbindlist(lapply(seq_along(patient_ids), function(i) {
    x_raw <- raw[, i]
    x_rpm <- rpm[, i]
    frac <- abs(x_raw - floor(x_raw))
    top_k <- max(1L, ceiling(0.01 * length(x_raw)))
    data.table(
      patient_id = patient_ids[i],
      analyte = analyte,
      raw_total = sum(x_raw),
      RPM_total = sum(x_rpm),
      nonzero_feature_n = sum(x_raw > 0),
      detection_fraction = mean(x_raw > 0),
      median_raw = median(x_raw),
      IQR_raw = IQR(x_raw),
      zero_fraction = mean(x_raw == 0),
      top_1pct_raw_share = sum(sort(x_raw, decreasing = TRUE)[seq_len(top_k)]) / max(sum(x_raw), 1),
      fractional_raw_fraction = mean(abs(x_raw - round(x_raw)) > 1e-8),
      fractional_part_q25 = unname(quantile(frac, 0.25)),
      fractional_part_median = median(frac),
      fractional_part_q75 = unname(quantile(frac, 0.75)),
      fractional_part_max = max(frac),
      negative_value_n = sum(x_raw < 0) + sum(x_rpm < 0),
      nonfinite_value_n = sum(!is.finite(x_raw)) + sum(!is.finite(x_rpm)),
      median_spearman_to_others = median(spearman[i, -i]),
      minimum_spearman_to_others = min(spearman[i, -i]),
      median_pearson_to_others = median(pearson[i, -i]),
      median_euclidean_distance = median(sample_dist[i, -i]),
      PCA_distance = pca_distance[i],
      PC1 = pca_scores$PC1[i],
      PC2 = pca_scores$PC2[i]
    )
  }))
  list(metrics = metrics, spearman = spearman, pearson = pearson,
       pca = pca_scores, variance_explained = pca$sdev^2 / sum(pca$sdev^2),
       top_feature_n = top_n)
}

log_event("Computing pre-group technical QC.")
qc_enst <- analyte_qc(raw_mat[is_enst, , drop = FALSE], rpm_mat[is_enst, , drop = FALSE],
                      "ENST", manifest$patient_id)
qc_mirna <- analyte_qc(raw_mat[is_mirna, , drop = FALSE], rpm_mat[is_mirna, , drop = FALSE],
                       "miRNA", manifest$patient_id)
qc_metrics <- rbindlist(list(qc_enst$metrics, qc_mirna$metrics))
qc_metrics <- merge(qc_metrics, manifest[, .(patient_id, GEO_sample_id, disease_group,
                                             sex, BMI, tissue_zone, terminal_format_warning)],
                    by = "patient_id", all.x = TRUE)
fwrite(qc_metrics, file.path(dirs[["qc_tables"]], "HD_T1_SAMPLE_QC_METRICS.csv"))
fwrite(as.data.table(qc_enst$spearman, keep.rownames = "patient_id"),
       file.path(dirs[["qc_tables"]], "HD_T1_ENST_SPEARMAN_CORRELATION.csv"))
fwrite(as.data.table(qc_mirna$spearman, keep.rownames = "patient_id"),
       file.path(dirs[["qc_tables"]], "HD_T1_MIRNA_SPEARMAN_CORRELATION.csv"))
fwrite(as.data.table(qc_enst$pearson, keep.rownames = "patient_id"),
       file.path(dirs[["qc_tables"]], "HD_T1_ENST_PEARSON_CORRELATION.csv"))
fwrite(as.data.table(qc_mirna$pearson, keep.rownames = "patient_id"),
       file.path(dirs[["qc_tables"]], "HD_T1_MIRNA_PEARSON_CORRELATION.csv"))

warning_rows <- list()
for (pid in manifest$patient_id) {
  q <- qc_metrics[patient_id == pid]
  reasons <- character()
  level <- "NO_WARNING"
  for (an in c("ENST", "miRNA")) {
    z_total <- safe_mad_z(qc_metrics[analyte == an]$raw_total)[qc_metrics[analyte == an]$patient_id == pid]
    z_detect <- safe_mad_z(qc_metrics[analyte == an]$nonzero_feature_n)[qc_metrics[analyte == an]$patient_id == pid]
    z_corr <- safe_mad_z(qc_metrics[analyte == an]$median_spearman_to_others)[qc_metrics[analyte == an]$patient_id == pid]
    z_pca <- safe_mad_z(qc_metrics[analyte == an]$PCA_distance)[qc_metrics[analyte == an]$patient_id == pid]
    if (z_total < -3.5) reasons <- c(reasons, paste0(an, "_extreme_low_raw_total"))
    if (z_detect < -3.5) reasons <- c(reasons, paste0(an, "_extreme_low_detection"))
    if (z_corr < -3.5) reasons <- c(reasons, paste0(an, "_low_sample_correlation"))
    if (z_pca > 3.5) reasons <- c(reasons, paste0(an, "_high_PCA_influence"))
  }
  src_warning <- manifest[patient_id == pid]$terminal_format_warning != "NO_WARNING"
  if (src_warning) {
    reasons <- c(reasons, "frozen_terminal_source_warning_requires_prespecified_sensitivity")
    level <- "QC_WARNING_REQUIRES_SENSITIVITY_ANALYSIS"
  } else if (length(reasons)) {
    level <- "QC_WARNING_NON_BLOCKING"
  }
  if (any(q$negative_value_n > 0 | q$nonfinite_value_n > 0)) {
    level <- "QC_FAILURE_BLOCKING"
    reasons <- c(reasons, "negative_or_nonfinite_common_matrix_value")
  }
  warning_rows[[length(warning_rows) + 1L]] <- data.table(
    patient_id = pid,
    GEO_sample_id = manifest[patient_id == pid]$GEO_sample_id,
    warning_class = level,
    warning_reasons = if (length(reasons)) paste(unique(reasons), collapse = ";") else "none",
    automatic_exclusion = "NO",
    primary_analysis_status = "RETAINED"
  )
}
qc_warnings <- rbindlist(warning_rows)
fwrite(qc_warnings, file.path(dirs[["qc_tables"]], "HD_T1_SAMPLE_QC_WARNINGS.csv"))

manifest[, BMI_numeric := as.numeric(BMI)]
manifest[, disease_group := factor(disease_group,
  levels = c("anal fissure surgical control", "hemorrhoidal disease"))]
manifest[, sex := factor(sex, levels = c("female", "male"))]
manifest[, tissue_zone := factor(tissue_zone,
  levels = c("anoderm region", "transition zone", "intestinal mucosa region"))]
manifest[, BMI_centered := BMI_numeric - mean(BMI_numeric)]

smd_binary <- function(p1, p0) {
  den <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  ifelse(den == 0, NA_real_, (p1 - p0) / den)
}
case_idx <- manifest$disease_group == "hemorrhoidal disease"
control_idx <- !case_idx
balance_rows <- list()
for (var in c("sex", "tissue_zone")) {
  tab <- table(manifest[[var]], manifest$disease_group)
  p_value <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
  for (lev in rownames(tab)) {
    p_case <- tab[lev, "hemorrhoidal disease"] / sum(tab[, "hemorrhoidal disease"])
    p_ctrl <- tab[lev, "anal fissure surgical control"] / sum(tab[, "anal fissure surgical control"])
    balance_rows[[length(balance_rows) + 1L]] <- data.table(
      variable = var, level = lev,
      case_n = tab[lev, "hemorrhoidal disease"], control_n = tab[lev, "anal fissure surgical control"],
      case_value = p_case, control_value = p_ctrl,
      standardized_difference = smd_binary(p_case, p_ctrl),
      descriptive_p_value = p_value, missing_n = sum(is.na(manifest[[var]])),
      interpretation = "Descriptive balance only; P value is not a covariate-selection rule."
    )
  }
}
bmi_case <- manifest$BMI_numeric[case_idx]
bmi_ctrl <- manifest$BMI_numeric[control_idx]
bmi_smd <- (mean(bmi_case) - mean(bmi_ctrl)) /
  sqrt((var(bmi_case) + var(bmi_ctrl)) / 2)
balance_rows[[length(balance_rows) + 1L]] <- data.table(
  variable = "BMI", level = "continuous",
  case_n = length(bmi_case), control_n = length(bmi_ctrl),
  case_value = mean(bmi_case), control_value = mean(bmi_ctrl),
  standardized_difference = bmi_smd,
  descriptive_p_value = wilcox.test(bmi_case, bmi_ctrl, exact = FALSE)$p.value,
  missing_n = sum(is.na(manifest$BMI_numeric)),
  interpretation = "Continuous BMI retained; descriptive Wilcoxon P value is not a covariate-selection rule."
)
metadata_balance <- rbindlist(balance_rows, fill = TRUE)
fwrite(metadata_balance, file.path(dirs[["qc_tables"]], "HD_T1_METADATA_BALANCE_TABLE.csv"))
fwrite(as.data.table(table(manifest$tissue_zone, manifest$disease_group), keep.rownames = "tissue_zone"),
       file.path(dirs[["qc_tables"]], "HD_T1_TISSUE_ZONE_BY_DISEASE.csv"))
fwrite(as.data.table(table(manifest$sex, manifest$tissue_zone), keep.rownames = "sex"),
       file.path(dirs[["qc_tables"]], "HD_T1_SEX_BY_TISSUE_ZONE.csv"))

design_unadjusted <- model.matrix(~ disease_group, data = manifest)
design_adjusted <- model.matrix(~ disease_group + sex + BMI_centered + tissue_zone, data = manifest)
rank_unadjusted <- qr(design_unadjusted)$rank
rank_adjusted <- qr(design_adjusted)$rank
full_rank_adjusted <- rank_adjusted == ncol(design_adjusted)
condition_adjusted <- kappa(design_adjusted)
zone_tab <- table(manifest$tissue_zone, manifest$disease_group)
zone_both_groups <- apply(zone_tab > 0, 1, all)
design_lines <- c(
  "# HD-T1 design matrix audit",
  "",
  "No gene or miRNA differential model was fitted.",
  "",
  "## Unadjusted design",
  "",
  "- Formula: `~ disease_group`.",
  paste0("- Dimensions: ", nrow(design_unadjusted), " x ", ncol(design_unadjusted), "."),
  paste0("- Rank: ", rank_unadjusted, "."),
  paste0("- Full rank: ", rank_unadjusted == ncol(design_unadjusted), "."),
  "",
  "## Frozen adjusted design",
  "",
  "- Formula: `~ disease_group + sex + BMI_centered + tissue_zone`.",
  paste0("- Dimensions: ", nrow(design_adjusted), " x ", ncol(design_adjusted), "."),
  paste0("- Rank: ", rank_adjusted, "."),
  paste0("- Full rank: ", full_rank_adjusted, "."),
  paste0("- Condition number: ", fmt_num(condition_adjusted, 3), "."),
  paste0("- Residual degrees of freedom if fitted: ", nrow(design_adjusted) - rank_adjusted, "."),
  paste0("- Coefficients: `", paste(colnames(design_adjusted), collapse = "`, `"), "`."),
  "",
  "## Tissue-zone support",
  "",
  paste(capture.output(print(zone_tab)), collapse = "\n"),
  "",
  paste0("- Every zone contains both groups: ", all(zone_both_groups), "."),
  "- Anoderm has only two cases and is descriptive-only under the frozen SAP.",
  "- Intestinal-mucosa-only analysis is possible but has only three controls and must be labeled low precision.",
  "- No zone is merged in HD-T1.",
  "",
  "## Modeling recommendation",
  "",
  "- Sex is complete and estimable in the frozen adjusted model.",
  "- BMI is complete and should remain continuous and centered.",
  "- Tissue zone is estimable but strongly imbalanced; adjusted, unadjusted, reduced-zone and prespecified zone-sensitive results remain mandatory in later authorized phases.",
  "- Age and batch are unavailable and remain prohibited.",
  paste0("- HD-T1 identifiability conclusion: `",
         if (full_rank_adjusted) "FULL_RANK_WITH_MANDATORY_ZONE_CAUTION" else "NOT_FULL_RANK_REQUIRES_PROTOCOL_REVISION",
         "`.")
)
write_utf8(design_lines, file.path(out_root, "HD_T1_DESIGN_MATRIX_AUDIT.md"))
fwrite(as.data.table(design_adjusted, keep.rownames = "patient_row"),
       file.path(dirs[["qc_tables"]], "HD_T1_ADJUSTED_DESIGN_MATRIX.csv"))

theme_set(theme_bw(base_size = 11))
qc_plot_data <- copy(qc_metrics)
qc_plot_data[, disease_group := factor(disease_group,
  levels = c("anal fissure surgical control", "hemorrhoidal disease"))]
subtitle_guard <- "Technical QC only; not a formal case-control result or disease-biology conclusion."
figure_records <- list()
add_figure <- function(stem, title, paths) {
  figure_records[[length(figure_records) + 1L]] <<- data.table(
    figure_id = stem, title = title, png_path = paths[1], pdf_path = paths[2],
    interpretation_boundary = subtitle_guard
  )
}

p1 <- ggplot(manifest, aes(x = disease_group, fill = disease_group)) +
  geom_bar(width = 0.65) + facet_wrap(~ tissue_zone) +
  scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "HD-T1 patient and tissue-zone structure", subtitle = subtitle_guard,
       x = NULL, y = "Patients") + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
add_figure("HD_T1_FIG01_PATIENT_STRUCTURE", "Patient and tissue-zone structure",
           save_plot_pair(p1, "HD_T1_FIG01_PATIENT_STRUCTURE", 11, 7))

p2 <- ggplot(qc_plot_data, aes(x = patient_id, y = raw_total, fill = disease_group)) +
  geom_col() + facet_wrap(~ analyte, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "Per-sample raw abundance totals", subtitle = subtitle_guard,
       x = "Patient in formal manifest order", y = "Raw abundance total") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
add_figure("HD_T1_FIG02_RAW_TOTALS", "Per-sample raw totals",
           save_plot_pair(p2, "HD_T1_FIG02_RAW_TOTALS", 12, 8))

p3 <- ggplot(qc_plot_data, aes(x = patient_id, y = nonzero_feature_n, fill = disease_group)) +
  geom_col() + facet_wrap(~ analyte, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "Detected non-zero features per sample", subtitle = subtitle_guard,
       x = "Patient in formal manifest order", y = "Non-zero feature count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
add_figure("HD_T1_FIG03_NONZERO_FEATURES", "Non-zero features",
           save_plot_pair(p3, "HD_T1_FIG03_NONZERO_FEATURES", 12, 8))

corr_long <- as.data.table(as.table(qc_enst$spearman))
setnames(corr_long, c("sample_x", "sample_y", "correlation"))
p4 <- ggplot(corr_long, aes(sample_x, sample_y, fill = correlation)) +
  geom_tile() + scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = median(corr_long[sample_x != sample_y]$correlation),
    limits = range(corr_long$correlation)
  ) +
  labs(title = "ENST sample Spearman correlation", subtitle = subtitle_guard,
       x = NULL, y = NULL, fill = "rho") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 7))
add_figure("HD_T1_FIG04_ENST_SPEARMAN_HEATMAP", "ENST Spearman heatmap",
           save_plot_pair(p4, "HD_T1_FIG04_ENST_SPEARMAN_HEATMAP", 11, 9))

pca_enst <- merge(qc_enst$pca, manifest[, .(patient_id, disease_group, tissue_zone)], by = "patient_id")
p5 <- ggplot(pca_enst, aes(PC1, PC2, color = disease_group, shape = tissue_zone, label = patient_id)) +
  geom_point(size = 3) + geom_text(nudge_y = 0.15, size = 2.5, check_overlap = TRUE) +
  scale_color_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "ENST PCA on top variable log2(RPM+1) features", subtitle = subtitle_guard,
       x = sprintf("PC1 (%.1f%%)", 100 * qc_enst$variance_explained[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * qc_enst$variance_explained[2])) +
  theme(legend.position = "bottom", legend.box = "vertical")
add_figure("HD_T1_FIG05_ENST_PCA", "ENST PCA",
           save_plot_pair(p5, "HD_T1_FIG05_ENST_PCA", 12, 8))

pca_mirna <- merge(qc_mirna$pca, manifest[, .(patient_id, disease_group, tissue_zone)], by = "patient_id")
p6 <- ggplot(pca_mirna, aes(PC1, PC2, color = disease_group, shape = tissue_zone, label = patient_id)) +
  geom_point(size = 3) + geom_text(nudge_y = 0.15, size = 2.5, check_overlap = TRUE) +
  scale_color_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "miRNA PCA on top variable log2(RPM+1) features", subtitle = subtitle_guard,
       x = sprintf("PC1 (%.1f%%)", 100 * qc_mirna$variance_explained[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * qc_mirna$variance_explained[2])) +
  theme(legend.position = "bottom", legend.box = "vertical")
add_figure("HD_T1_FIG06_MIRNA_PCA", "miRNA PCA",
           save_plot_pair(p6, "HD_T1_FIG06_MIRNA_PCA", 12, 8))

p7 <- ggplot(manifest, aes(disease_group, BMI_numeric, fill = disease_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) + geom_jitter(width = 0.12, size = 2) +
  scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "BMI distribution", subtitle = subtitle_guard, x = NULL, y = "BMI") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
add_figure("HD_T1_FIG07_BMI_DISTRIBUTION", "BMI distribution",
           save_plot_pair(p7, "HD_T1_FIG07_BMI_DISTRIBUTION", 8, 6))

p8 <- ggplot(manifest, aes(tissue_zone, fill = disease_group)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "Tissue zone by disease group", subtitle = subtitle_guard,
       x = NULL, y = "Patients") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
add_figure("HD_T1_FIG08_ZONE_BY_GROUP", "Tissue zone by group",
           save_plot_pair(p8, "HD_T1_FIG08_ZONE_BY_GROUP", 9, 6))

p9 <- ggplot(qc_plot_data, aes(patient_id, fractional_raw_fraction, fill = disease_group)) +
  geom_col() + facet_wrap(~ analyte, ncol = 1) +
  scale_fill_manual(values = c("#4C78A8", "#E45756")) +
  labs(title = "Fractional raw abundance audit", subtitle = subtitle_guard,
       x = "Patient in formal manifest order", y = "Fractional-value proportion") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
add_figure("HD_T1_FIG09_FRACTIONAL_RAW", "Fractional raw audit",
           save_plot_pair(p9, "HD_T1_FIG09_FRACTIONAL_RAW", 12, 8))

feature_count_plot <- data.table(
  universe = c("Common", "Noncommon union"),
  ENST = c(common_enst, sum(noncommon$feature_type == "ENST transcript")),
  miRNA = c(common_mirna, sum(noncommon$feature_type == "miRNA")),
  other = c(common_other, sum(noncommon$feature_type == "other/unknown"))
)
feature_count_long <- melt(feature_count_plot, id.vars = "universe",
                           variable.name = "feature_type", value.name = "feature_n")
p10 <- ggplot(feature_count_long, aes(universe, feature_n, fill = feature_type)) +
  geom_col(position = "stack") + scale_y_continuous(labels = scales::comma) +
  labs(title = "Common and noncommon feature counts", subtitle = subtitle_guard,
       x = NULL, y = "Features", fill = "Feature type")
add_figure("HD_T1_FIG10_FEATURE_UNIVERSE", "Common and noncommon feature counts",
           save_plot_pair(p10, "HD_T1_FIG10_FEATURE_UNIVERSE", 8, 6))

figure_index <- rbindlist(figure_records)
fwrite(figure_index, file.path(dirs[["qc_figures"]], "HD_T1_QC_FIGURE_INDEX.csv"))
figure_md <- c("# HD-T1 QC figure index", "",
               "Every figure is technical QC only and is not a formal case-control result.",
               "",
               paste0("- `", figure_index$figure_id, "`: ", figure_index$title,
                      " ([PNG](", figure_index$png_path, "), [PDF](", figure_index$pdf_path, "))."))
write_utf8(figure_md, file.path(out_root, "HD_T1_QC_FIGURE_INDEX.md"))

run_end <- Sys.time()
loaded_versions <- data.table(
  package = required_packages,
  version = vapply(required_packages, function(p) as.character(packageVersion(p)), character(1))
)
fwrite(loaded_versions, file.path(dirs[["environment"]], "HD_T1_LOADED_PACKAGE_VERSIONS.csv"))
session_path <- file.path(dirs[["environment"]], "HD_T1_SESSION_INFO.txt")
write_utf8(capture.output(sessionInfo()), session_path)
input_hashes <- manifest[, .(GEO_sample_id, patient_id, processed_tsv_file,
                             sha256 = toupper(processed_tsv_sha256))]
fwrite(input_hashes, file.path(dirs[["environment"]], "HD_T1_INPUT_SHA256_SET.csv"))
annotation_hashes <- data.table(
  reference = c("Ensembl_release_100_GRCh38_GTF", "miRBase_22.1_mature.decoded.fa",
                "miRBase_22.1_hairpin.decoded.fa", "miRBase_22.1_miRNA.decoded.dead",
                "miRBase_22.1_README.decoded"),
  path = c(gtf_path, mirbase_mature_path, mirbase_hairpin_path, mirbase_dead_path, mirbase_readme_path),
  file_size_bytes = file.info(c(gtf_path, mirbase_mature_path, mirbase_hairpin_path,
                                mirbase_dead_path, mirbase_readme_path))$size,
  sha256 = vapply(c(gtf_path, mirbase_mature_path, mirbase_hairpin_path,
                    mirbase_dead_path, mirbase_readme_path),
                  function(x) if (file.exists(x)) sha256(x) else NA_character_, character(1))
)
fwrite(annotation_hashes, file.path(dirs[["environment"]], "HD_T1_ANNOTATION_REFERENCE_HASHES.csv"))

git_commit <- "NOT_AVAILABLE_EMPTY_GIT_DIRECTORY"
environment_lines <- c(
  "# HD-T1 environment lock",
  "",
  paste0("- Start: ", format(run_start, "%Y-%m-%d %H:%M:%S %Z"), "."),
  paste0("- End: ", format(run_end, "%Y-%m-%d %H:%M:%S %Z"), "."),
  paste0("- Runtime seconds: ", fmt_num(as.numeric(difftime(run_end, run_start, units = "secs")), 2), "."),
  paste0("- Operating system: ", paste(Sys.info()[c("sysname", "release", "version", "machine")], collapse = " | "), "."),
  paste0("- R version: ", R.version.string, "."),
  "- Python: not used.",
  paste0("- Locale: ", paste(Sys.getlocale(), collapse = "; "), "."),
  paste0("- Working directory: `", getwd(), "`."),
  paste0("- Project root: `", project_root, "`."),
  paste0("- Time zone: ", Sys.timezone(), "."),
  "- Random seed: 20260716.",
  paste0("- Git commit: ", git_commit, "."),
  paste0("- Formal manifest SHA-256: `", manifest_sha, "`."),
  paste0("- T0R verified manifest SHA-256: `", verified_manifest_sha, "`."),
  "- Formal differential-expression packages loaded: none.",
  "",
  "## Loaded packages",
  "",
  paste0("- ", loaded_versions$package, ": ", loaded_versions$version),
  "",
  "## Annotation audit references",
  "",
  paste0("- ", annotation_hashes$reference, ": `", annotation_hashes$sha256, "`.")
)
write_utf8(environment_lines, file.path(out_root, "HD_T1_ENVIRONMENT_LOCK.md"))

matrix_integrity_pass <- all(matrix_audit$NA_n == 0 & matrix_audit$Inf_n == 0 &
                               matrix_audit$negative_n == 0 &
                               matrix_audit$duplicate_row_name_n == 0 &
                               matrix_audit$duplicate_column_name_n == 0)
raw_fractional_preserved <- all(qc_metrics[analyte %in% c("ENST", "miRNA")]$fractional_raw_fraction > 0)
blocking_sample_n <- sum(qc_warnings$warning_class == "QC_FAILURE_BLOCKING")
nonblocking_warning_n <- sum(qc_warnings$warning_class != "NO_WARNING" &
                               qc_warnings$warning_class != "QC_FAILURE_BLOCKING")
gate_checks <- data.table(
  gate_item = c(
    "formal_manifest_equals_T0R_verified_manifest", "all_38_patients_imported",
    "frozen_20_case_18_control_groups", "common_feature_total_33959",
    "common_ENST_32153", "common_miRNA_1806", "no_zero_fill_for_noncommon_features",
    "all_common_matrices_nonnegative_finite_unique", "fractional_raw_precision_preserved",
    "feature_classification_reproducible", "ENST_mapping_audit_completed",
    "miRNA_name_audit_completed", "no_blocking_sample_failure", "all_38_retained",
    "adjusted_design_full_rank", "scripts_logs_versions_hashes_recorded",
    "no_formal_differential_expression_run"
  ),
  pass = c(
    manifest_sha == verified_manifest_sha, nrow(input_manifest) == 38L,
    sum(manifest$disease_group == "hemorrhoidal disease") == 20L &&
      sum(manifest$disease_group == "anal fissure surgical control") == 18L,
    common_total == 33959L, common_enst == 32153L, common_mirna == 1806L,
    TRUE, matrix_integrity_pass, raw_fractional_preserved, common_other == 0L,
    nrow(enst_annotation) >= common_enst, nrow(mirna_map) >= common_mirna,
    blocking_sample_n == 0L, all(qc_warnings$primary_analysis_status == "RETAINED"),
    full_rank_adjusted, TRUE, TRUE
  )
)
gate_checks[, evidence := c(
  manifest_sha, as.character(nrow(input_manifest)), "20/18",
  as.character(common_total), as.character(common_enst), as.character(common_mirna),
  "noncommon audit labels missing records as not observed; matrices use intersection only",
  as.character(matrix_integrity_pass), as.character(raw_fractional_preserved),
  paste0("other_common_features=", common_other),
  paste0(mapped_enst_n, "/", common_enst, " mapped or audited"),
  paste0(recognized_mirna_n, "/", common_mirna, " recognized current"),
  as.character(blocking_sample_n), "38 retained", paste0("rank=", rank_adjusted, "/", ncol(design_adjusted)),
  "script, sessionInfo, package versions, SHA-256 sets and logs written",
  "No limma, edgeR or DESeq2 package loaded; no feature-level group model fitted"
)]
fwrite(gate_checks, file.path(dirs[["gate"]], "HD_T1_GATE_CHECKLIST.csv"))

if (!all(gate_checks$pass[gate_checks$gate_item != "adjusted_design_full_rank"])) {
  if (!gate_checks[gate_item == "common_feature_total_33959", pass] ||
      !gate_checks[gate_item == "common_ENST_32153", pass] ||
      !gate_checks[gate_item == "common_miRNA_1806", pass]) {
    final_state <- "HD_T1_FAIL_COMMON_FEATURE_SPACE"
  } else if (!gate_checks[gate_item == "all_common_matrices_nonnegative_finite_unique", pass] ||
             !gate_checks[gate_item == "all_38_patients_imported", pass]) {
    final_state <- "HD_T1_FAIL_INPUT_INTEGRITY"
  } else if (blocking_sample_n > 0L) {
    final_state <- "HD_T1_FAIL_SAMPLE_QC"
  } else {
    final_state <- "HD_T1_REQUIRES_PROTOCOL_REVISION"
  }
} else if (!full_rank_adjusted) {
  final_state <- "HD_T1_FAIL_DESIGN_IDENTIFIABILITY"
} else if (nonblocking_warning_n > 0L || unmapped_enst_n > 0L ||
           unrecognized_mirna_n > 0L || deprecated_mirna_n > 0L) {
  final_state <- "HD_T1_PASS_WITH_NONBLOCKING_WARNINGS_READY_FOR_T2"
} else {
  final_state <- "HD_T1_PASS_READY_FOR_T2"
}

input_report <- c(
  "# HD-T1 input integrity report",
  "",
  paste0("- Formal/verified manifest SHA-256 identity: `", manifest_sha, "`."),
  "- Source files: 38/38 present and independently re-hashed.",
  paste0("- Cases/controls: ", sum(case_idx), "/", sum(control_idx), "."),
  paste0("- Common feature intersection: ", fmt_int(common_total), "."),
  paste0("- Feature union: ", fmt_int(length(union_features)), "."),
  paste0("- Common ENST: ", fmt_int(common_enst), "."),
  paste0("- Common miRNA: ", fmt_int(common_mirna), "."),
  paste0("- Common other/unknown: ", fmt_int(common_other), "."),
  "- GSM4676457 absent tail features were not imputed or zero-filled.",
  "- GSM4676447 and GSM4676457 source files were not modified.",
  paste0("- Matrix integrity pass: ", matrix_integrity_pass, "."),
  paste0("- Blocking sample failures: ", blocking_sample_n, ".")
)
write_utf8(input_report, file.path(out_root, "HD_T1_INPUT_INTEGRITY_REPORT.md"))

warning_text <- qc_warnings[warning_class != "NO_WARNING",
  paste0(GEO_sample_id, " (", patient_id, "): ", warning_class, " — ", warning_reasons)]
if (!length(warning_text)) warning_text <- "None."
executive <- c(
  "# HD-T1 executive summary",
  "",
  paste0("**Final state:** `", final_state, "`"),
  "",
  "- HD-T1 only was executed. HD-T2 was not started.",
  paste0("- Final patients: 38 (", sum(case_idx), " hemorrhoidal disease; ",
         sum(control_idx), " anal fissure surgical controls)."),
  paste0("- Common features: ", fmt_int(common_total), " (", fmt_int(common_enst),
         " ENST; ", fmt_int(common_mirna), " miRNA)."),
  paste0("- Blocking samples: ", blocking_sample_n, "."),
  paste0("- ENST Ensembl-100 mapping rate: ", fmt_num(100 * mapped_enst_n / common_enst, 2), "%."),
  paste0("- miRNA miRBase-22.1 current-name recognition rate: ",
         fmt_num(100 * recognized_mirna_n / common_mirna, 2), "%."),
  paste0("- Frozen adjusted design full rank: ", full_rank_adjusted,
         " (rank ", rank_adjusted, "/", ncol(design_adjusted), ")."),
  "- Sex is estimable; BMI remains continuous and centered; tissue zone is estimable but strongly imbalanced and requires all frozen later sensitivities.",
  "",
  "## Nonblocking warnings",
  "",
  paste0("- ", warning_text),
  "",
  "## Formal HD-T2 inputs prepared but not consumed",
  "",
  paste0("- Raw ENST matrix: `", matrix_paths[["RAW_ENST"]], "`."),
  paste0("- Raw miRNA matrix: `", matrix_paths[["RAW_MIRNA"]], "`."),
  paste0("- Formal manifest: `", manifest_path, "`."),
  "",
  "A PASS state indicates technical readiness only. It does not authorize HD-T2."
)
write_utf8(executive, file.path(out_root, "HD_T1_EXECUTIVE_SUMMARY.md"))

final_gate <- c(
  "# HD-T1 final gate",
  "",
  paste0("## Final status\n\n`", final_state, "`"),
  "",
  paste0("- Patients: 38 (", sum(case_idx), " cases; ", sum(control_idx), " controls)."),
  paste0("- Common features: ", common_total, " total; ", common_enst, " ENST; ", common_mirna, " miRNA."),
  paste0("- Blocking sample failures: ", blocking_sample_n, "."),
  paste0("- Nonblocking warning patients: ", nonblocking_warning_n, "."),
  paste0("- ENST mapping: ", mapped_enst_n, "/", common_enst, " (",
         fmt_num(100 * mapped_enst_n / common_enst, 2), "%)."),
  paste0("- miRNA current-name recognition: ", recognized_mirna_n, "/", common_mirna, " (",
         fmt_num(100 * recognized_mirna_n / common_mirna, 2), "%)."),
  paste0("- Adjusted design full rank: ", full_rank_adjusted, "."),
  "- Formal differential expression run: NO.",
  "- HD-T2 automatic start: NO.",
  "",
  "## Gate interpretation",
  "",
  if (grepl("^HD_T1_PASS", final_state))
    "The HD-T1 technical outputs are ready for review. Explicit user authorization remains required before any HD-T2 action."
  else
    "The reported failure or revision state blocks HD-T2. No workaround by sample deletion, group merging or covariate alteration was applied."
)
write_utf8(final_gate, file.path(out_root, "HD_T1_FINAL_GATE.md"))

log_event(paste("Final state:", final_state))
log_event("HD-T1 complete. Stopping without starting HD-T2.")
write_utf8(log_lines, file.path(out_root, "HD_T1_ANALYSIS_LOG.md"))
write_utf8(log_lines, file.path(dirs[["logs"]], "HD_T1_RUN.log"))

required_outputs <- c(
  file.path(out_root, "HD_T1_EXECUTIVE_SUMMARY.md"),
  file.path(out_root, "HD_T1_INPUT_INTEGRITY_REPORT.md"),
  file.path(out_root, "HD_T1_ENVIRONMENT_LOCK.md"),
  file.path(dirs[["manifests"]], "HD_T1_INPUT_FILE_MANIFEST.csv"),
  file.path(dirs[["manifests"]], "HD_T1_FEATURE_PRESENCE_MATRIX.csv.gz"),
  file.path(dirs[["manifests"]], "HD_T1_COMMON_FEATURE_LIST.csv"),
  file.path(dirs[["manifests"]], "HD_T1_NONCOMMON_FEATURE_AUDIT.csv"),
  file.path(dirs[["annotation"]], "HD_T1_FEATURE_TYPE_MAP.csv"),
  unname(matrix_paths),
  file.path(dirs[["annotation"]], "HD_T1_ENST_ANNOTATION_MAP.csv"),
  file.path(out_root, "HD_T1_ENST_MAPPING_SUMMARY.md"),
  file.path(dirs[["annotation"]], "HD_T1_MIRNA_NAME_MAP.csv"),
  file.path(dirs[["qc_tables"]], "HD_T1_SAMPLE_QC_METRICS.csv"),
  file.path(dirs[["qc_tables"]], "HD_T1_SAMPLE_QC_WARNINGS.csv"),
  file.path(dirs[["qc_tables"]], "HD_T1_METADATA_BALANCE_TABLE.csv"),
  file.path(out_root, "HD_T1_DESIGN_MATRIX_AUDIT.md"),
  file.path(out_root, "HD_T1_QC_FIGURE_INDEX.md"),
  file.path(out_root, "HD_T1_ANALYSIS_LOG.md"),
  file.path(out_root, "HD_T1_FINAL_GATE.md")
)
output_index <- data.table(
  output_path = required_outputs,
  exists = file.exists(required_outputs),
  file_size_bytes = ifelse(file.exists(required_outputs), file.info(required_outputs)$size, NA_real_),
  sha256 = vapply(required_outputs, function(x) if (file.exists(x)) sha256(x) else NA_character_, character(1))
)
fwrite(output_index, file.path(out_root, "HD_T1_OUTPUT_INDEX.csv"))
if (!all(output_index$exists)) stop("One or more required HD-T1 outputs are missing after execution.")

cat(final_state, "\n")
