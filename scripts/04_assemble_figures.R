options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(scales)
  library(svglite)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "results", "HD_T7B")

dirs <- c(
  "00_logs", "01_source_provenance", "02_figure_data", "03_main_figures",
  "04_main_tables", "05_supplementary_figures", "06_supplementary_tables",
  "07_legends", "08_visual_qa", "09_environment", "10_gate"
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (d in dirs) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out, "00_logs", "HD_T7B_ASSEMBLY.log")
zz <- file(log_file, open = "wt", encoding = "UTF-8")
sink(zz, type = "output", split = TRUE)
sink(zz, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("HD-T7B assembly started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
cat("Root:", root, "\n")
cat("Backend: R only\n")

rel <- function(...) file.path(...)
abs_path <- function(...) file.path(root, ...)
out_path <- function(...) file.path(out, ...)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_text_utf8 <- function(lines, path) {
  writeLines(enc2utf8(lines), con = path, useBytes = TRUE)
}

sha256_file <- function(path) {
  toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

assert_equal <- function(observed, expected, label, tolerance = 0) {
  ok <- if (is.numeric(observed) && is.numeric(expected)) {
    isTRUE(all.equal(as.numeric(observed), as.numeric(expected), tolerance = tolerance))
  } else {
    identical(as.character(observed), as.character(expected))
  }
  if (!ok) stop(sprintf("NUMERIC INCONSISTENCY: %s observed=%s expected=%s", label, paste(observed, collapse = "|"), paste(expected, collapse = "|")))
  invisible(TRUE)
}

required_packages <- c("ggplot2", "gridExtra", "scales", "svglite", "data.table", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

palette <- c(
  case = "#4C78A8",
  control = "#E39C37",
  unadjusted = "#697386",
  no_zone = "#4C78A8",
  full = "#B65A65",
  reduced = "#7B6BA8",
  attenuated = "#3D9C89",
  reversal = "#D65F5F",
  neutral = "#B9C0C8",
  dark = "#2B2B2B",
  light = "#EEF1F4",
  low_precision = "#6E9DC5",
  sample_sensitive = "#D6A04B",
  zone_dependent = "#B56B8C",
  model_dependent = "#8A7DB8",
  unstable = "#8C8C8C",
  only_unadjusted = "#C9C9C9"
)

model_levels <- c(
  "Completely\nunadjusted",
  "Sex/BMI adjusted,\nno zone",
  "Full adjusted,\nincluding zone",
  "Reduced-zone"
)
model_colors <- setNames(
  c(palette[["unadjusted"]], palette[["no_zone"]], palette[["full"]], palette[["reduced"]]),
  model_levels
)

theme_hd <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(3.2, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#444444"),
      plot.caption = element_text(size = base_size - 0.8, colour = "#555555", hjust = 0),
      plot.margin = margin(4, 5, 4, 5),
      panel.grid = element_blank()
    )
}
theme_set(theme_hd())

as_grob <- function(x) {
  if (inherits(x, "ggplot")) ggplotGrob(x) else x
}

tag_grob <- function(x, tag) {
  grobTree(
    as_grob(x),
    textGrob(tag, x = unit(1.5, "mm"), y = unit(1, "npc") - unit(1.5, "mm"),
             just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = 8, fontface = "bold"))
  )
}

draw_grob <- function(x) {
  grid.newpage()
  grid.draw(as_grob(x))
}

save_grob_bundle <- function(x, stem, width_mm, height_mm, dpi = 600) {
  pdf_file <- paste0(stem, ".pdf")
  svg_file <- paste0(stem, ".svg")
  png_file <- paste0(stem, "_600DPI.png")
  grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", onefile = TRUE)
  draw_grob(x)
  dev.off()
  svglite::svglite(svg_file, width = width_mm / 25.4, height = height_mm / 25.4, system_fonts = list(sans = "Arial"))
  draw_grob(x)
  dev.off()
  grDevices::png(
    png_file, width = width_mm, height = height_mm, units = "mm", res = dpi,
    type = "cairo", family = "Arial", bg = "white"
  )
  draw_grob(x)
  dev.off()
  invisible(c(pdf_file, svg_file, png_file))
}

save_panel <- function(x, fig, panel, width_mm = 89, height_mm = 70) {
  stem <- out_path("03_main_figures", sprintf("HD_T7B_FIGURE_%d_PANEL_%s", fig, panel))
  save_grob_bundle(tag_grob(x, tolower(panel)), stem, width_mm, height_mm)
}

copy_canonical <- function(source, destination_name) {
  dest <- out_path(destination_name)
  if (!file.copy(source, dest, overwrite = TRUE)) stop("Failed to copy canonical output: ", source)
  dest
}

box_plot <- function(boxes, arrows = NULL, notes = NULL, xlim = c(0, 10), ylim = c(0, 10)) {
  p <- ggplot() +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(5, 7, 5, 7))
  if (!is.null(arrows) && nrow(arrows)) {
    p <- p + geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.45, colour = "#6E7781",
      arrow = arrow(length = unit(2.2, "mm"), type = "closed")
    )
  }
  p <- p +
    geom_label(
      data = boxes,
      aes(x = x, y = y, label = label, fill = fill),
      family = "Arial", size = 2.45, label.size = 0.3, label.r = unit(1.5, "mm"),
      colour = palette[["dark"]], lineheight = 0.95, label.padding = unit(1.8, "mm")
    ) +
    scale_fill_identity()
  if (!is.null(notes) && nrow(notes)) {
    p <- p + geom_text(
      data = notes, aes(x = x, y = y, label = label),
      family = "Arial", size = 2.25, colour = "#4F5660", lineheight = 0.95
    )
  }
  p
}

label_count <- function(x) scales::comma(x, accuracy = 1)
prop_pct <- function(x, denom) 100 * x / denom

# -------------------------------------------------------------------------
# Frozen inputs
# -------------------------------------------------------------------------

files <- list(
  manifest = abs_path("HD_DATASET_AND_SAMPLE_MANIFEST.csv"),
  r3_coverage = abs_path("HD_R3_ORIGINAL_PAPER_COVERAGE_MATRIX.csv"),
  r4_claims = abs_path("HD_R4_FINAL_CLAIM_EVIDENCE_MATRIX.csv"),
  t1_balance = abs_path("results", "HD_T1", "04_qc_tables", "HD_T1_METADATA_BALANCE_TABLE.csv"),
  t1_zone = abs_path("results", "HD_T1", "04_qc_tables", "HD_T1_TISSUE_ZONE_BY_DISEASE.csv"),
  t1_warnings = abs_path("results", "HD_T1", "04_qc_tables", "HD_T1_SAMPLE_QC_WARNINGS.csv"),
  t1_design = abs_path("results", "HD_T1", "04_qc_tables", "HD_T1_ADJUSTED_DESIGN_MATRIX.csv"),
  t2_unadjusted = abs_path("results", "HD_T2", "HD_T2_GENE_DE_UNADJUSTED.csv.gz"),
  t2_no_zone = abs_path("results", "HD_T2", "HD_T2_GENE_DE_NO_ZONE_ADJUSTMENT.csv.gz"),
  t2_adjusted = abs_path("results", "HD_T2", "HD_T2_GENE_DE_ADJUSTED.csv.gz"),
  t2_reduced = abs_path("results", "HD_T2", "HD_T2_GENE_DE_REDUCED_ZONE.csv.gz"),
  t2_stability = abs_path("results", "HD_T2", "HD_T2_GENE_STABILITY_CLASSIFICATION.csv"),
  t2_loo = abs_path("results", "HD_T2", "HD_T2_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv"),
  t2_influence = abs_path("results", "HD_T2", "HD_T2_SAMPLE_INFLUENCE_METRICS.csv"),
  t2_ledger = abs_path("results", "HD_T2", "HD_T2_RESULTS_LEDGER.md.csv"),
  t2_modules = abs_path("results", "HD_T2", "HD_T2_MODULE_STABILITY_RESULTS.csv"),
  t2_module_scores = abs_path("results", "HD_T2", "HD_T2_MODULE_DIFFERENTIAL_RESULTS.csv"),
  t3_filter = abs_path("results", "HD_T3", "02_filtering", "HD_T3_FILTER_COUNTS_BY_FAMILY.csv"),
  t3_unadjusted = abs_path("results", "HD_T3", "HD_T3_MIRNA_DE_UNADJUSTED.csv.gz"),
  t3_no_zone = abs_path("results", "HD_T3", "HD_T3_MIRNA_DE_NO_ZONE_ADJUSTMENT.csv.gz"),
  t3_adjusted = abs_path("results", "HD_T3", "HD_T3_MIRNA_DE_ADJUSTED.csv.gz"),
  t3_reduced = abs_path("results", "HD_T3", "HD_T3_MIRNA_DE_REDUCED_ZONE.csv.gz"),
  t3_stability = abs_path("results", "HD_T3", "HD_T3_MIRNA_STABILITY_CLASSIFICATION.csv"),
  t3_zone_audit = abs_path("results", "HD_T3", "HD_T3_ZONE_DEPENDENCE_AUDIT.csv"),
  t3_loo = abs_path("results", "HD_T3", "HD_T3_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv"),
  t3_influence = abs_path("results", "HD_T3", "HD_T3_SAMPLE_INFLUENCE_METRICS.csv"),
  t3_ledger = abs_path("results", "HD_T3", "06_tables", "HD_T3_RESULTS_LEDGER.csv")
)

missing_inputs <- names(files)[!file.exists(unlist(files))]
if (length(missing_inputs)) stop("Missing frozen inputs: ", paste(missing_inputs, collapse = ", "))

manifest <- read_csv(files$manifest)
t1_balance <- read_csv(files$t1_balance)
t1_zone <- read_csv(files$t1_zone)
t1_warnings <- read_csv(files$t1_warnings)
t1_design <- read_csv(files$t1_design)
t2_unadj <- read_csv(files$t2_unadjusted)
t2_noz <- read_csv(files$t2_no_zone)
t2_adj <- read_csv(files$t2_adjusted)
t2_red <- read_csv(files$t2_reduced)
t2_stab <- read_csv(files$t2_stability)
t2_loo <- read_csv(files$t2_loo)
t2_infl <- read_csv(files$t2_influence)
t2_ledger <- read_csv(files$t2_ledger)
t2_modules <- read_csv(files$t2_modules)
t2_module_scores <- read_csv(files$t2_module_scores)
t3_filter <- read_csv(files$t3_filter)
t3_unadj <- read_csv(files$t3_unadjusted)
t3_noz <- read_csv(files$t3_no_zone)
t3_adj <- read_csv(files$t3_adjusted)
t3_red <- read_csv(files$t3_reduced)
t3_stab <- read_csv(files$t3_stability)
t3_zone <- read_csv(files$t3_zone_audit)
t3_loo <- read_csv(files$t3_loo)
t3_infl <- read_csv(files$t3_influence)
t3_ledger <- read_csv(files$t3_ledger)
r4_claims <- read_csv(files$r4_claims)

# -------------------------------------------------------------------------
# Numeric lock and data transformations
# -------------------------------------------------------------------------

case_label <- "hemorrhoidal disease"
control_label <- "anal fissure surgical control"
manifest$disease_group <- trimws(manifest$disease_group)

n_total <- nrow(manifest)
n_case <- sum(manifest$disease_group == case_label)
n_control <- sum(manifest$disease_group == control_label)
assert_equal(n_total, 38, "patients")
assert_equal(n_case, 20, "case patients")
assert_equal(n_control, 18, "control patients")
assert_equal(sum(manifest$mRNA_present == "YES"), 38, "mRNA complete")
assert_equal(sum(manifest$miRNA_present == "YES"), 38, "miRNA complete")

zone_counts <- as.data.frame(with(manifest, table(tissue_zone, disease_group)), stringsAsFactors = FALSE)
names(zone_counts) <- c("tissue_zone", "disease_group", "n")
zone_counts <- zone_counts[zone_counts$n > 0, ]
get_zone_n <- function(zone, group) zone_counts$n[zone_counts$tissue_zone == zone & zone_counts$disease_group == group]
assert_equal(get_zone_n("anoderm region", case_label), 2, "case anoderm")
assert_equal(get_zone_n("anoderm region", control_label), 10, "control anoderm")
assert_equal(get_zone_n("transition zone", case_label), 4, "case transition")
assert_equal(get_zone_n("transition zone", control_label), 5, "control transition")
assert_equal(get_zone_n("intestinal mucosa region", case_label), 14, "case intestinal")
assert_equal(get_zone_n("intestinal mucosa region", control_label), 3, "control intestinal")

gene_counts <- c(
  unadjusted = sum(t2_unadj$BH_FDR < 0.05, na.rm = TRUE),
  no_zone = sum(t2_noz$BH_FDR < 0.05, na.rm = TRUE),
  full = sum(t2_adj$BH_FDR < 0.05, na.rm = TRUE),
  reduced = sum(t2_red$BH_FDR < 0.05, na.rm = TRUE)
)
assert_equal(nrow(t2_adj), 12161, "tested genes")
assert_equal(unname(gene_counts), c(6111, 5389, 0, 0), "gene FDR ladder")

mirna_mature_count <- function(d) sum(d$feature_family == "MATURE", na.rm = TRUE)
mirna_fdr_count <- function(d) sum(d$feature_family == "MATURE" & d$BH_FDR < 0.05, na.rm = TRUE)
mirna_counts <- c(
  unadjusted = mirna_fdr_count(t3_unadj),
  no_zone = mirna_fdr_count(t3_noz),
  full = mirna_fdr_count(t3_adj),
  reduced = mirna_fdr_count(t3_red)
)
assert_equal(nrow(t3_adj), 697, "all filtered miRNA records")
assert_equal(mirna_mature_count(t3_adj), 670, "filtered mature miRNA")
assert_equal(unname(mirna_counts), c(46, 45, 0, 0), "mature miRNA FDR ladder")

gene_migration <- merge(
  t2_noz[, c("ensembl_gene_id", "gene_symbol", "log2_fold_change", "BH_FDR")],
  t2_adj[, c("ensembl_gene_id", "log2_fold_change", "BH_FDR")],
  by = "ensembl_gene_id", suffixes = c("_no_zone", "_adjusted"), all = FALSE
)
gene_migration$no_zone_significant <- gene_migration$BH_FDR_no_zone < 0.05
gene_migration$attenuated <- abs(gene_migration$log2_fold_change_adjusted) < abs(gene_migration$log2_fold_change_no_zone)
gene_migration$direction_reversed <- sign(gene_migration$log2_fold_change_adjusted) != sign(gene_migration$log2_fold_change_no_zone)
gene_sig <- gene_migration[gene_migration$no_zone_significant, ]
gene_attenuated <- sum(gene_sig$attenuated)
gene_reversed <- sum(gene_sig$direction_reversed)
gene_before_median <- median(abs(gene_sig$log2_fold_change_no_zone))
gene_after_median <- median(abs(gene_sig$log2_fold_change_adjusted))
assert_equal(gene_attenuated, 5255, "attenuated no-zone significant genes")
assert_equal(gene_reversed, 1016, "gene direction reversals")
assert_equal(round(gene_before_median, 3), 0.513, "gene median absolute log2FC before")
assert_equal(round(gene_after_median, 3), 0.149, "gene median absolute log2FC after")

mirna_noz_m <- t3_noz[t3_noz$feature_family == "MATURE", ]
mirna_adj_m <- t3_adj[t3_adj$feature_family == "MATURE", ]
mirna_migration <- merge(
  mirna_noz_m[, c("original_miRNA_id", "normalized_miRNA_id", "logFC", "BH_FDR")],
  mirna_adj_m[, c("original_miRNA_id", "logFC", "BH_FDR")],
  by = "original_miRNA_id", suffixes = c("_no_zone", "_adjusted"), all = FALSE
)
mirna_migration$no_zone_significant <- mirna_migration$BH_FDR_no_zone < 0.05
mirna_migration$attenuated <- abs(mirna_migration$logFC_adjusted) < abs(mirna_migration$logFC_no_zone)
mirna_migration$direction_reversed <- sign(mirna_migration$logFC_adjusted) != sign(mirna_migration$logFC_no_zone)
mirna_sig <- mirna_migration[mirna_migration$no_zone_significant, ]
mirna_attenuated <- sum(mirna_sig$attenuated)
mirna_reversed <- sum(mirna_sig$direction_reversed)
mirna_before_median <- median(abs(mirna_sig$logFC_no_zone))
mirna_after_median <- median(abs(mirna_sig$logFC_adjusted))
assert_equal(mirna_attenuated, 41, "attenuated no-zone significant mature miRNAs")
assert_equal(mirna_reversed, 5, "mature miRNA direction reversals")
assert_equal(round(mirna_before_median, 3), 1.767, "miRNA median absolute log2FC before")
assert_equal(round(mirna_after_median, 3), 0.537, "miRNA median absolute log2FC after")

gene_class_counts <- as.data.frame(table(t2_stab$stability_class), stringsAsFactors = FALSE)
names(gene_class_counts) <- c("stability_class", "n")
gene_class_counts <- gene_class_counts[gene_class_counts$n > 0, ]
gene_class_map <- c(
  DIRECTIONALLY_STABLE_LOW_PRECISION = "Low-precision\ndirectional",
  SAMPLE_SENSITIVE = "Sample-sensitive",
  ZONE_DEPENDENT = "Zone-dependent",
  MODEL_DEPENDENT = "Model-dependent"
)
gene_class_counts$display_class <- unname(gene_class_map[gene_class_counts$stability_class])
gene_class_counts$analyte <- "mRNA genes"
gene_class_counts$denominator <- 12161
gene_class_counts$percent <- 100 * gene_class_counts$n / gene_class_counts$denominator
assert_equal(sum(t2_stab$stability_class == "DIRECTIONALLY_STABLE_LOW_PRECISION"), 4210, "low-precision genes")
assert_equal(sum(t2_stab$stability_class == "ROBUST"), 0, "robust genes")

t3_stab_filtered <- t3_stab[!is.na(t3_stab$primary_logFC), ]
t3_stab_mature <- t3_stab_filtered[t3_stab_filtered$feature_family == "MATURE", ]
mirna_class_counts <- as.data.frame(table(t3_stab_mature$stability_class), stringsAsFactors = FALSE)
names(mirna_class_counts) <- c("stability_class", "n")
mirna_class_counts <- mirna_class_counts[mirna_class_counts$n > 0, ]
mirna_class_map <- c(
  DIRECTIONALLY_STABLE_LOW_PRECISION = "Low-precision\ndirectional",
  SAMPLE_SENSITIVE = "Sample-sensitive",
  ZONE_DEPENDENT = "Zone-dependent",
  MODEL_DEPENDENT = "Model-dependent",
  UNSTABLE = "Unstable",
  ONLY_UNADJUSTED = "Only unadjusted"
)
mirna_class_counts$display_class <- unname(mirna_class_map[mirna_class_counts$stability_class])
mirna_class_counts$analyte <- "Mature miRNAs"
mirna_class_counts$denominator <- 670
mirna_class_counts$percent <- 100 * mirna_class_counts$n / mirna_class_counts$denominator
assert_equal(sum(t3_stab_mature$stability_class == "DIRECTIONALLY_STABLE_LOW_PRECISION"), 229, "low-precision mature miRNAs")
assert_equal(sum(t3_stab_mature$stability_class == "SAMPLE_SENSITIVE"), 333, "sample-sensitive mature miRNAs")
assert_equal(sum(t3_stab_mature$stability_class == "ZONE_DEPENDENT"), 93, "zone-dependent mature miRNAs")
assert_equal(sum(t3_stab_mature$stability_class == "MODEL_DEPENDENT"), 10, "model-dependent mature miRNAs")
assert_equal(sum(t3_stab_mature$stability_class == "UNSTABLE"), 4, "unstable mature miRNAs")
assert_equal(sum(t3_stab_mature$stability_class == "ONLY_UNADJUSTED"), 1, "only-unadjusted mature miRNAs")
assert_equal(sum(t3_stab_mature$adjusted_stable), 0, "adjusted-stable mature miRNAs")
assert_equal(sum(t3_stab_filtered$stability_class == "SAMPLE_SENSITIVE"), 348, "all-family sample-sensitive miRNAs")
assert_equal(sum(t3_stab_filtered$stability_class == "ZONE_DEPENDENT"), 95, "all-family zone-dependent miRNAs")

stable_camera <- sum(t2_modules$stable_predefined_module)
score_fdr <- sum(t2_modules$score_BH_FDR < 0.05, na.rm = TRUE)
assert_equal(stable_camera, 3, "stable camera modules")
assert_equal(score_fdr, 0, "FDR-significant patient-level module scores")

source_inventory <- data.frame(
  source_key = names(files),
  source_file = vapply(files, function(x) substring(normalizePath(x, winslash = "/", mustWork = TRUE), nchar(root) + 2L), character(1)),
  source_file_sha256 = vapply(files, sha256_file, character(1)),
  readable = "YES",
  frozen_input = "YES",
  stringsAsFactors = FALSE
)
write_csv_utf8(source_inventory, out_path("01_source_provenance", "HD_T7B_SOURCE_FILE_INVENTORY.csv"))

# -------------------------------------------------------------------------
# Figure source-data files
# -------------------------------------------------------------------------

fig1_flow <- data.frame(
  item = c("Patients", "Hemorrhoidal disease", "Anal fissure surgical controls", "Paired mRNA", "Paired miRNA", "Genes tested", "All miRNA records after filtering", "Mature miRNAs tested"),
  value = c(38, 20, 18, 38, 38, 12161, 697, 670),
  source = c("manifest", "manifest", "manifest", "manifest", "manifest", "T2 adjusted results", "T3 filter ledger", "T3 mature results")
)
write_csv_utf8(fig1_flow, out_path("02_figure_data", "HD_T7B_FIGURE_1A_SOURCE_DATA.csv"))
write_csv_utf8(zone_counts, out_path("02_figure_data", "HD_T7B_FIGURE_1C_SOURCE_DATA.csv"))

model_ladder <- data.frame(
  model = factor(model_levels, levels = model_levels),
  formula = c(
    "expression ~ disease group",
    "expression ~ disease group + sex + centered BMI",
    "expression ~ disease group + sex + centered BMI + expression-derived zone",
    "exclude anoderm; retain sex, centered BMI and zone"
  ),
  n_patients = c(38, 38, 38, 26)
)
write_csv_utf8(model_ladder, out_path("02_figure_data", "HD_T7B_FIGURE_1D_SOURCE_DATA.csv"))

gene_ladder <- data.frame(
  model = factor(model_levels, levels = model_levels),
  fdr_significant_n = unname(gene_counts),
  analyte = "mRNA genes"
)
write_csv_utf8(gene_ladder, out_path("02_figure_data", "HD_T7B_FIGURE_2A_SOURCE_DATA.csv"))
write_csv_utf8(gene_migration, out_path("02_figure_data", "HD_T7B_FIGURE_2B_SOURCE_DATA.csv"))
gene_migration_summary <- data.frame(
  metric = c("No-zone FDR signals", "Attenuated after zone addition", "Direction reversals", "Median |log2FC| before", "Median |log2FC| after"),
  value = c(5389, gene_attenuated, gene_reversed, gene_before_median, gene_after_median),
  denominator = c(12161, 5389, 5389, 5389, 5389)
)
write_csv_utf8(gene_migration_summary, out_path("02_figure_data", "HD_T7B_FIGURE_2C_SOURCE_DATA.csv"))
write_csv_utf8(gene_class_counts, out_path("02_figure_data", "HD_T7B_FIGURE_2D_STABILITY_SOURCE_DATA.csv"))
write_csv_utf8(t2_infl, out_path("02_figure_data", "HD_T7B_FIGURE_2D_INFLUENCE_SOURCE_DATA.csv"))

mirna_filter_long <- data.frame(
  family = factor(c("Mature", "Precursor", "Unresolved"), levels = c("Mature", "Precursor", "Unresolved")),
  before = c(1684, 92, 30),
  after = c(670, 14, 13)
)
write_csv_utf8(mirna_filter_long, out_path("02_figure_data", "HD_T7B_FIGURE_3A_FILTER_SOURCE_DATA.csv"))
mirna_ladder <- data.frame(
  model = factor(model_levels, levels = model_levels),
  fdr_significant_n = unname(mirna_counts),
  analyte = "Mature miRNAs"
)
write_csv_utf8(mirna_ladder, out_path("02_figure_data", "HD_T7B_FIGURE_3A_MODEL_SOURCE_DATA.csv"))
write_csv_utf8(mirna_migration, out_path("02_figure_data", "HD_T7B_FIGURE_3B_SOURCE_DATA.csv"))
mirna_migration_summary <- data.frame(
  metric = c("No-zone FDR signals", "Attenuated after zone addition", "Direction reversals", "Median |log2FC| before", "Median |log2FC| after"),
  value = c(45, mirna_attenuated, mirna_reversed, mirna_before_median, mirna_after_median),
  denominator = c(670, 45, 45, 45, 45)
)
write_csv_utf8(mirna_migration_summary, out_path("02_figure_data", "HD_T7B_FIGURE_3C_SOURCE_DATA.csv"))
write_csv_utf8(mirna_class_counts, out_path("02_figure_data", "HD_T7B_FIGURE_3D_STABILITY_SOURCE_DATA.csv"))
write_csv_utf8(t3_infl, out_path("02_figure_data", "HD_T7B_FIGURE_3D_INFLUENCE_SOURCE_DATA.csv"))

cross_ladder <- rbind(
  transform(gene_ladder, analyte = "mRNA genes"),
  transform(mirna_ladder, analyte = "Mature miRNAs")
)
write_csv_utf8(cross_ladder, out_path("02_figure_data", "HD_T7B_FIGURE_4A_SOURCE_DATA.csv"))
cross_migration <- data.frame(
  analyte = rep(c("mRNA genes", "Mature miRNAs"), each = 2),
  metric = rep(c("Attenuated", "Direction reversed"), 2),
  numerator = c(gene_attenuated, gene_reversed, mirna_attenuated, mirna_reversed),
  denominator = c(5389, 5389, 45, 45)
)
cross_migration$percent <- 100 * cross_migration$numerator / cross_migration$denominator
write_csv_utf8(cross_migration, out_path("02_figure_data", "HD_T7B_FIGURE_4B_SOURCE_DATA.csv"))
cross_classes <- rbind(gene_class_counts, mirna_class_counts)
write_csv_utf8(cross_classes, out_path("02_figure_data", "HD_T7B_FIGURE_4C_SOURCE_DATA.csv"))

# -------------------------------------------------------------------------
# Figure 1
# -------------------------------------------------------------------------

fig1a_boxes <- data.frame(
  x = c(2.0, 4.7, 7.8, 7.8, 7.8),
  y = c(7.2, 7.2, 8.5, 6.2, 2.8),
  label = c(
    "38 patients\n20 disease | 18 controls",
    "One patient-level\nCombo-Seq library",
    "mRNA branch\n12,161 genes tested",
    "miRNA branch\n697 retained records\n670 mature tested",
    "Frozen four-model ladder\nand sensitivity audit"
  ),
  fill = c("#E8EEF6", "#F1F3F5", "#DCE8F4", "#F6E8D7", "#EAE5F2")
)
fig1a_arrows <- data.frame(
  x = c(2.8, 5.8, 5.8, 7.8, 7.8),
  y = c(7.2, 7.2, 7.2, 7.7, 5.4),
  xend = c(3.6, 6.7, 6.7, 7.8, 7.8),
  yend = c(7.2, 8.25, 6.45, 4.0, 4.0)
)
p1a <- box_plot(fig1a_boxes, fig1a_arrows) +
  ggtitle("Patient-level paired-analyte design")

fig1b_boxes <- data.frame(
  x = c(1.8, 5.0, 8.2, 5.0),
  y = c(7.5, 7.5, 7.5, 3.2),
  label = c(
    "Keratinocyte/sebocyte\nexpression signatures",
    "ssGSEA + clustering\n(original study)",
    "Deposited three-level\nzone classification",
    "Expression-derived covariate\nnot independent pathology or anatomy"
  ),
  fill = c("#E8EEF6", "#EAE5F2", "#F6E8D7", "#FBE3E3")
)
fig1b_arrows <- data.frame(
  x = c(3.0, 6.2, 8.2),
  y = c(7.5, 7.5, 6.4),
  xend = c(3.8, 7.0, 5.9),
  yend = c(7.5, 7.5, 3.8)
)
p1b <- box_plot(fig1b_boxes, fig1b_arrows) +
  ggtitle("Provenance of the zone classification")

zone_plot <- zone_counts
zone_plot$tissue_zone <- factor(
  zone_plot$tissue_zone,
  levels = c("anoderm region", "transition zone", "intestinal mucosa region"),
  labels = c("Anoderm", "Transition", "Intestinal mucosa")
)
zone_plot$disease_group <- factor(
  zone_plot$disease_group,
  levels = c(case_label, control_label),
  labels = c("Hemorrhoidal disease", "Anal fissure controls")
)
p1c <- ggplot(zone_plot, aes(tissue_zone, n, fill = disease_group)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72), vjust = -0.35, size = 2.5, family = "Arial") +
  scale_fill_manual(values = c("Hemorrhoidal disease" = palette[["case"]], "Anal fissure controls" = palette[["control"]])) +
  scale_y_continuous(breaks = seq(0, 16, 4), limits = c(0, 16), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "Sparse and asymmetric disease-by-zone overlap",
    x = "Expression-derived histological zone", y = "Patients",
    fill = NULL,
    caption = "Zone was inferred from expression signatures, ssGSEA and clustering."
  ) +
  theme_hd() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "top")

fig1d_boxes <- data.frame(
  x = rep(5, 4),
  y = c(8.3, 6.2, 4.1, 1.8),
  label = c(
    "1  Completely unadjusted\nexpression ~ disease group  |  n = 38",
    "2  Sex/BMI adjusted, no zone\n+ sex + centered BMI  |  n = 38",
    "3  Full adjusted, including zone\n+ expression-derived zone  |  n = 38",
    "4  Reduced-zone sensitivity\nexclude anoderm; retain zone  |  n = 26"
  ),
  fill = c("#E4E7EC", "#D8E7F4", "#F4DDE0", "#E7E0F2")
)
fig1d_arrows <- data.frame(
  x = c(5, 5, 5),
  y = c(7.45, 5.35, 3.15),
  xend = c(5, 5, 5),
  yend = c(6.95, 4.85, 2.65)
)
p1d <- box_plot(fig1d_boxes, fig1d_arrows) +
  ggtitle("Four distinct model contexts")

save_panel(p1a, 1, "A")
save_panel(p1b, 1, "B")
save_panel(p1c, 1, "C")
save_panel(p1d, 1, "D")

fig1 <- arrangeGrob(
  tag_grob(p1a, "a"), tag_grob(p1b, "b"),
  tag_grob(p1c, "c"), tag_grob(p1d, "d"),
  ncol = 2, widths = c(1, 1), heights = c(1, 1)
)
fig1_stem <- out_path("03_main_figures", "HD_T7B_FIGURE_1")
save_grob_bundle(fig1, fig1_stem, 183, 142)
copy_canonical(paste0(fig1_stem, ".pdf"), "HD_T7B_FIGURE_1.pdf")
copy_canonical(paste0(fig1_stem, "_600DPI.png"), "HD_T7B_FIGURE_1_600DPI.png")
copy_canonical(paste0(fig1_stem, ".svg"), "HD_T7B_FIGURE_1.svg")

# -------------------------------------------------------------------------
# Figure 2
# -------------------------------------------------------------------------

p2a <- ggplot(gene_ladder, aes(model, fdr_significant_n, fill = model)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = comma(fdr_significant_n)), vjust = ifelse(gene_ladder$fdr_significant_n == 0, -0.6, -0.3), size = 2.7, family = "Arial") +
  scale_fill_manual(values = model_colors, guide = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 6800), breaks = c(0, 2000, 4000, 6000), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "FDR-significant genes across model contexts",
    subtitle = "Primary zone-addition comparison: 5,389 → 0",
    x = NULL, y = "Genes at BH FDR < 0.05",
    caption = "6,111 is completely unadjusted context, not the isolated zone-addition contrast."
  ) +
  theme_hd() +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

gene_migration$plot_class <- ifelse(
  gene_migration$no_zone_significant & gene_migration$direction_reversed, "Direction reversed",
  ifelse(gene_migration$no_zone_significant, "No-zone FDR signal", "Other tested gene")
)
gene_plot_order <- order(gene_migration$plot_class != "Other tested gene")
gene_scatter <- gene_migration[gene_plot_order, ]
p2b <- ggplot(gene_scatter, aes(log2_fold_change_no_zone, log2_fold_change_adjusted, colour = plot_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#C8CDD2") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#C8CDD2") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "#6E7781") +
  geom_point(size = 0.42, alpha = 0.38, stroke = 0) +
  scale_colour_manual(values = c(
    "Other tested gene" = "#C9CED4",
    "No-zone FDR signal" = palette[["no_zone"]],
    "Direction reversed" = palette[["reversal"]]
  )) +
  coord_equal() +
  labs(
    title = "Gene-level disease effects contract after zone addition",
    x = "Sex/BMI-adjusted no-zone log2FC",
    y = "Full adjusted log2FC",
    colour = NULL
  ) +
  theme_hd() +
  theme(legend.position = "top")

gene_change_plot <- data.frame(
  metric = factor(c("Attenuated", "Direction reversed"), levels = c("Attenuated", "Direction reversed")),
  n = c(gene_attenuated, gene_reversed),
  percent = c(prop_pct(gene_attenuated, 5389), prop_pct(gene_reversed, 5389))
)
p2c <- ggplot(gene_change_plot, aes(metric, percent, fill = metric)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%s / 5,389\n(%.1f%%)", comma(n), percent)), vjust = -0.25, size = 2.65, family = "Arial", lineheight = 0.9) +
  scale_fill_manual(values = c("Attenuated" = palette[["attenuated"]], "Direction reversed" = palette[["reversal"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Migration among no-zone FDR signals",
    subtitle = sprintf("Median |log2FC| %.3f → %.3f", gene_before_median, gene_after_median),
    x = NULL, y = "Proportion of 5,389 no-zone signals"
  ) +
  theme_hd()

gene_class_counts$display_class <- factor(
  gene_class_counts$display_class,
  levels = c("Low-precision\ndirectional", "Sample-sensitive", "Zone-dependent", "Model-dependent")
)
p2d_top <- ggplot(gene_class_counts, aes(display_class, n, fill = display_class)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = comma(n)), vjust = -0.25, size = 2.35, family = "Arial") +
  scale_fill_manual(values = c(
    "Low-precision\ndirectional" = palette[["low_precision"]],
    "Sample-sensitive" = palette[["sample_sensitive"]],
    "Zone-dependent" = palette[["zone_dependent"]],
    "Model-dependent" = palette[["model_dependent"]]
  ), guide = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 4700), breaks = c(0, 2000, 4000), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Frozen mutually exclusive stability classes", x = NULL, y = "Genes") +
  theme_hd(base_size = 6.6) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

t2_infl$patient_order <- rank(t2_infl$maximum_absolute_log2FC_change, ties.method = "first")
p2d_bottom <- ggplot(t2_infl, aes(patient_order, maximum_absolute_log2FC_change)) +
  geom_hline(yintercept = median(t2_infl$maximum_absolute_log2FC_change), linetype = 2, linewidth = 0.35, colour = "#6E7781") +
  geom_point(aes(colour = warning_class), size = 1.15, alpha = 0.9) +
  scale_colour_manual(values = c(
    "NO_WARNING" = "#7E8791",
    "QC_WARNING_NON_BLOCKING" = "#D6A04B",
    "QC_WARNING_REQUIRES_SENSITIVITY_ANALYSIS" = "#B65A65"
  )) +
  labs(
    x = "Patients ranked by leave-one-out effect change", y = "Maximum |Δlog2FC|", colour = NULL,
    caption = "Grey: no warning; gold: non-blocking warning; rose: sensitivity-required warning."
  ) +
  theme_hd(base_size = 6.2) +
  theme(legend.position = "none")

p2d <- arrangeGrob(p2d_top, p2d_bottom, ncol = 1, heights = c(1.2, 0.85))

save_panel(p2a, 2, "A")
save_panel(p2b, 2, "B")
save_panel(p2c, 2, "C")
save_panel(p2d, 2, "D")

fig2 <- arrangeGrob(
  tag_grob(p2a, "a"), tag_grob(p2b, "b"),
  tag_grob(p2c, "c"), tag_grob(p2d, "d"),
  ncol = 2, widths = c(0.9, 1.1), heights = c(1, 1.05)
)
fig2_stem <- out_path("03_main_figures", "HD_T7B_FIGURE_2")
save_grob_bundle(fig2, fig2_stem, 183, 150)
copy_canonical(paste0(fig2_stem, ".pdf"), "HD_T7B_FIGURE_2.pdf")
copy_canonical(paste0(fig2_stem, "_600DPI.png"), "HD_T7B_FIGURE_2_600DPI.png")
copy_canonical(paste0(fig2_stem, ".svg"), "HD_T7B_FIGURE_2.svg")

# -------------------------------------------------------------------------
# Figure 3
# -------------------------------------------------------------------------

mirna_filter_plot <- rbind(
  data.frame(family = mirna_filter_long$family, stage = "Before filtering", n = mirna_filter_long$before),
  data.frame(family = mirna_filter_long$family, stage = "After filtering", n = mirna_filter_long$after)
)
p3a_top <- ggplot(mirna_filter_plot, aes(family, n, fill = stage)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = comma(n)), position = position_dodge(width = 0.72), vjust = -0.25, size = 2.2, family = "Arial") +
  scale_fill_manual(values = c("Before filtering" = "#BFC7D2", "After filtering" = palette[["no_zone"]])) +
  scale_y_continuous(labels = comma, limits = c(0, 1900), breaks = c(0, 500, 1000, 1500), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "1,806 deposited records → 697 retained", x = NULL, y = "miRNA records", fill = NULL) +
  theme_hd(base_size = 6.4) +
  theme(legend.position = "top")

p3a_bottom <- ggplot(mirna_ladder, aes(model, fdr_significant_n, fill = model)) +
  geom_col(width = 0.66, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = fdr_significant_n), vjust = ifelse(mirna_ladder$fdr_significant_n == 0, -0.6, -0.3), size = 2.3, family = "Arial") +
  scale_fill_manual(values = model_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 52), breaks = c(0, 15, 30, 45), expand = expansion(mult = c(0, 0.04))) +
  labs(subtitle = "Mature-miRNA FDR ladder: 46 → 45 → 0 → 0", x = NULL, y = "BH FDR < 0.05") +
  theme_hd(base_size = 6.2) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))
p3a <- arrangeGrob(p3a_top, p3a_bottom, ncol = 1, heights = c(1, 1))

mirna_migration$plot_class <- ifelse(
  mirna_migration$no_zone_significant & mirna_migration$direction_reversed, "Direction reversed",
  ifelse(mirna_migration$no_zone_significant, "No-zone FDR signal", "Other mature miRNA")
)
mirna_plot_order <- order(mirna_migration$plot_class != "Other mature miRNA")
mirna_scatter <- mirna_migration[mirna_plot_order, ]
p3b <- ggplot(mirna_scatter, aes(logFC_no_zone, logFC_adjusted, colour = plot_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#C8CDD2") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#C8CDD2") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "#6E7781") +
  geom_point(size = 0.75, alpha = 0.58, stroke = 0) +
  scale_colour_manual(values = c(
    "Other mature miRNA" = "#C9CED4",
    "No-zone FDR signal" = palette[["no_zone"]],
    "Direction reversed" = palette[["reversal"]]
  )) +
  coord_equal() +
  labs(
    title = "Mature-miRNA disease effects contract after zone addition",
    x = "Sex/BMI-adjusted no-zone log2FC",
    y = "Full adjusted log2FC",
    colour = NULL
  ) +
  theme_hd() +
  theme(legend.position = "top")

mirna_change_plot <- data.frame(
  metric = factor(c("Attenuated", "Direction reversed"), levels = c("Attenuated", "Direction reversed")),
  n = c(mirna_attenuated, mirna_reversed),
  percent = c(prop_pct(mirna_attenuated, 45), prop_pct(mirna_reversed, 45))
)
p3c <- ggplot(mirna_change_plot, aes(metric, percent, fill = metric)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%s / 45\n(%.1f%%)", comma(n), percent)), vjust = -0.25, size = 2.65, family = "Arial", lineheight = 0.9) +
  scale_fill_manual(values = c("Attenuated" = palette[["attenuated"]], "Direction reversed" = palette[["reversal"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Migration among no-zone mature-miRNA signals",
    subtitle = sprintf("Median |log2FC| %.3f → %.3f", mirna_before_median, mirna_after_median),
    x = NULL, y = "Proportion of 45 no-zone signals"
  ) +
  theme_hd()

mirna_class_counts$display_class <- factor(
  mirna_class_counts$display_class,
  levels = c("Low-precision\ndirectional", "Sample-sensitive", "Zone-dependent", "Model-dependent", "Unstable", "Only unadjusted")
)
p3d_top <- ggplot(mirna_class_counts, aes(display_class, n, fill = display_class)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.25, size = 2.2, family = "Arial") +
  scale_fill_manual(values = c(
    "Low-precision\ndirectional" = palette[["low_precision"]],
    "Sample-sensitive" = palette[["sample_sensitive"]],
    "Zone-dependent" = palette[["zone_dependent"]],
    "Model-dependent" = palette[["model_dependent"]],
    "Unstable" = palette[["unstable"]],
    "Only unadjusted" = palette[["only_unadjusted"]]
  ), guide = "none") +
  scale_y_continuous(limits = c(0, 370), breaks = c(0, 100, 200, 300), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Mature-miRNA stability classes (n = 670)", x = NULL, y = "Mature miRNAs") +
  theme_hd(base_size = 6.2) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

t3_infl$patient_order <- rank(t3_infl$maximum_absolute_logFC_change, ties.method = "first")
p3d_bottom <- ggplot(t3_infl, aes(patient_order, maximum_absolute_logFC_change)) +
  geom_hline(yintercept = median(t3_infl$maximum_absolute_logFC_change), linetype = 2, linewidth = 0.35, colour = "#6E7781") +
  geom_point(aes(colour = high_disease_effect_influence), size = 1.2, alpha = 0.9) +
  scale_colour_manual(values = c("FALSE" = "#7E8791", "TRUE" = palette[["reversal"]]), labels = c("FALSE" = "Not flagged", "TRUE" = "High influence")) +
  labs(x = "Patients ranked by leave-one-out effect change", y = "Maximum |Δlog2FC|", colour = NULL) +
  theme_hd(base_size = 6.2) +
  theme(legend.position = "top")
p3d <- arrangeGrob(p3d_top, p3d_bottom, ncol = 1, heights = c(1.2, 0.85))

save_panel(p3a, 3, "A")
save_panel(p3b, 3, "B")
save_panel(p3c, 3, "C")
save_panel(p3d, 3, "D")

fig3 <- arrangeGrob(
  tag_grob(p3a, "a"), tag_grob(p3b, "b"),
  tag_grob(p3c, "c"), tag_grob(p3d, "d"),
  ncol = 2, widths = c(0.92, 1.08), heights = c(1.05, 1)
)
fig3_stem <- out_path("03_main_figures", "HD_T7B_FIGURE_3")
save_grob_bundle(fig3, fig3_stem, 183, 154)
copy_canonical(paste0(fig3_stem, ".pdf"), "HD_T7B_FIGURE_3.pdf")
copy_canonical(paste0(fig3_stem, "_600DPI.png"), "HD_T7B_FIGURE_3_600DPI.png")
copy_canonical(paste0(fig3_stem, ".svg"), "HD_T7B_FIGURE_3.svg")

# -------------------------------------------------------------------------
# Figure 4
# -------------------------------------------------------------------------

cross_ladder$analyte <- factor(cross_ladder$analyte, levels = c("mRNA genes", "Mature miRNAs"))
p4a <- ggplot(cross_ladder, aes(model, fdr_significant_n, fill = model)) +
  geom_col(width = 0.66, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = comma(fdr_significant_n)), vjust = ifelse(cross_ladder$fdr_significant_n == 0, -0.55, -0.25), size = 2.2, family = "Arial") +
  facet_wrap(~analyte, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = model_colors, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Parallel significance migration",
    subtitle = "Independent y-scales; exact counts are printed",
    x = NULL, y = "Features at BH FDR < 0.05"
  ) +
  theme_hd(base_size = 6.3) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

cross_migration$metric <- factor(cross_migration$metric, levels = c("Attenuated", "Direction reversed"))
cross_migration$analyte <- factor(cross_migration$analyte, levels = c("mRNA genes", "Mature miRNAs"))
p4b <- ggplot(cross_migration, aes(analyte, percent, fill = metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", percent)), position = position_dodge(width = 0.72), vjust = -0.25, size = 2.3, family = "Arial") +
  scale_fill_manual(values = c("Attenuated" = palette[["attenuated"]], "Direction reversed" = palette[["reversal"]])) +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Effect migration among no-zone FDR signals",
    x = NULL, y = "Within-analyte proportion", fill = NULL,
    caption = "Denominators: 5,389 genes and 45 mature miRNAs."
  ) +
  theme_hd() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 10, hjust = 1))

cross_classes$display_class <- factor(
  cross_classes$display_class,
  levels = c("Low-precision\ndirectional", "Sample-sensitive", "Zone-dependent", "Model-dependent", "Unstable", "Only unadjusted")
)
cross_classes$analyte <- factor(cross_classes$analyte, levels = c("mRNA genes", "Mature miRNAs"))
p4c <- ggplot(cross_classes, aes(analyte, percent, fill = display_class)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(
    "Low-precision\ndirectional" = palette[["low_precision"]],
    "Sample-sensitive" = palette[["sample_sensitive"]],
    "Zone-dependent" = palette[["zone_dependent"]],
    "Model-dependent" = palette[["model_dependent"]],
    "Unstable" = palette[["unstable"]],
    "Only unadjusted" = palette[["only_unadjusted"]]
  ), drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(
    title = "Frozen evidence-class composition",
    x = NULL, y = "Percent of tested family", fill = "Mutually exclusive class",
    caption = "Rules were frozen separately for genes and miRNA families."
  ) +
  theme_hd(base_size = 6.4) +
  theme(legend.position = "right")

fig4d_boxes <- data.frame(
  x = rep(5, 5),
  y = c(8.6, 6.9, 5.2, 3.5, 1.4),
  label = c(
    "Apparent disease contrast\nmodels without expression-derived zone",
    "Add sex and centered BMI\nthen add deposited zone classification",
    "Full and reduced-zone FDR control\n0 genes; 0 mature miRNAs",
    "Residual directional evidence\nlow precision and patient sensitive",
    "Supported conclusion\ninferential sensitivity, not mechanism or biomarker"
  ),
  fill = c("#E8EEF6", "#F1F3F5", "#F6E8D7", "#EAE5F2", "#E4F0EA")
)
fig4d_arrows <- data.frame(
  x = rep(5, 4),
  y = c(7.95, 6.25, 4.55, 2.85),
  xend = rep(5, 4),
  yend = c(7.55, 5.85, 4.15, 2.25)
)
p4d <- box_plot(fig4d_boxes, fig4d_arrows) +
  ggtitle("Inference boundary")

save_panel(p4a, 4, "A")
save_panel(p4b, 4, "B")
save_panel(p4c, 4, "C")
save_panel(p4d, 4, "D")

fig4 <- arrangeGrob(
  tag_grob(p4a, "a"), tag_grob(p4b, "b"),
  tag_grob(p4c, "c"), tag_grob(p4d, "d"),
  ncol = 2, widths = c(1.08, 0.92), heights = c(1, 1)
)
fig4_stem <- out_path("03_main_figures", "HD_T7B_FIGURE_4")
save_grob_bundle(fig4, fig4_stem, 183, 148)
copy_canonical(paste0(fig4_stem, ".pdf"), "HD_T7B_FIGURE_4.pdf")
copy_canonical(paste0(fig4_stem, "_600DPI.png"), "HD_T7B_FIGURE_4_600DPI.png")
copy_canonical(paste0(fig4_stem, ".svg"), "HD_T7B_FIGURE_4.svg")

# -------------------------------------------------------------------------
# Main tables
# -------------------------------------------------------------------------

bmi_case <- manifest$BMI[manifest$disease_group == case_label]
bmi_control <- manifest$BMI[manifest$disease_group == control_label]
fmt_mean_sd <- function(x) sprintf("%.2f ± %.2f", mean(x), sd(x))

table1 <- data.frame(
  variable = c(
    "Patients", "Sex", "Sex", "BMI", "Expression-derived histological zone",
    "Expression-derived histological zone", "Expression-derived histological zone",
    "Paired analyte completeness", "Paired analyte completeness", "Age", "Batch",
    "Model hierarchy", "Model hierarchy", "Model hierarchy", "Model hierarchy"
  ),
  level_or_statistic = c(
    "Total n", "Female, n", "Male, n", "Mean ± SD",
    "Anoderm, n", "Transition, n", "Intestinal mucosa, n",
    "mRNA available, n", "miRNA available, n", "Not reported", "Not reported",
    "Completely unadjusted", "Sex/BMI adjusted, no zone",
    "Full adjusted including expression-derived zone", "Reduced-zone sensitivity"
  ),
  hemorrhoidal_disease = c(
    "20", "8", "12", fmt_mean_sd(bmi_case), "2", "4", "14", "20", "20",
    "Unavailable", "Unavailable", "n = 38", "n = 38", "n = 38", "n = 18 after excluding anoderm"
  ),
  anal_fissure_surgical_controls = c(
    "18", "9", "9", fmt_mean_sd(bmi_control), "10", "5", "3", "18", "18",
    "Unavailable", "Unavailable", "n = 38", "n = 38", "n = 38", "n = 8 after excluding anoderm"
  ),
  missingness = c(rep("0", 9), "38", "38", rep("Not applicable", 4)),
  interpretive_note = c(
    "Patient/library is the statistical unit.",
    "Descriptive; no P value used for covariate selection.",
    "Descriptive; no P value used for covariate selection.",
    "Continuous centered BMI was used in adjusted models.",
    "Expression-derived; not independently pathology confirmed.",
    "Expression-derived; not independently pathology confirmed.",
    "Expression-derived; not independently pathology confirmed.",
    "Same patient Combo-Seq library.",
    "Same patient Combo-Seq library.",
    "Not available in recovered public metadata.",
    "Not available in recovered public metadata.",
    "expression ~ disease group",
    "expression ~ disease group + sex + centered BMI",
    "expression ~ disease group + sex + centered BMI + expression-derived zone",
    "Anoderm excluded; transition and intestinal-mucosa classes retained with zone adjustment."
  )
)
table1_file <- out_path("04_main_tables", "HD_T7B_TABLE_1.csv")
write_csv_utf8(table1, table1_file)
copy_canonical(table1_file, "HD_T7B_TABLE_1.csv")

table2 <- data.frame(
  analyte_family = c("mRNA genes", "Mature miRNAs"),
  tested_features = c(12161, 670),
  completely_unadjusted_FDR_n = c(6111, 46),
  sex_BMI_adjusted_no_zone_FDR_n = c(5389, 45),
  full_adjusted_FDR_n = c(0, 0),
  reduced_zone_FDR_n = c(0, 0),
  no_zone_signals_attenuated_after_zone_n = c(gene_attenuated, mirna_attenuated),
  no_zone_signal_denominator = c(5389, 45),
  direction_reversal_n = c(gene_reversed, mirna_reversed),
  median_absolute_log2FC_before = c(gene_before_median, mirna_before_median),
  median_absolute_log2FC_after = c(gene_after_median, mirna_after_median),
  directionally_stable_low_precision_n = c(4210, 229),
  zone_or_model_dependent_n = c(
    sum(t2_stab$stability_class %in% c("ZONE_DEPENDENT", "MODEL_DEPENDENT")),
    sum(t3_stab_mature$stability_class %in% c("ZONE_DEPENDENT", "MODEL_DEPENDENT"))
  ),
  sample_sensitive_n = c(
    sum(t2_stab$stability_class == "SAMPLE_SENSITIVE"),
    sum(t3_stab_mature$stability_class == "SAMPLE_SENSITIVE")
  ),
  unstable_or_only_unadjusted_n = c(
    0,
    sum(t3_stab_mature$stability_class %in% c("UNSTABLE", "ONLY_UNADJUSTED"))
  ),
  robust_or_adjusted_stable_n = c(0, 0),
  principal_allowed_conclusion = c(
    "Gene-level disease inference is highly sensitive to adding the expression-derived zone classification; residual directional signals are low precision.",
    "Mature-miRNA inference shows parallel sensitivity without constituting regulatory integration or independent replication."
  ),
  principal_limitation = c(
    "Zone was inferred from mRNA expression; sparse within-zone overlap and patient influence limit disease-specific interpretation.",
    "Same cohort and expression-derived mRNA zone classification; mature, precursor and unresolved families remain separate."
  )
)
table2_file <- out_path("04_main_tables", "HD_T7B_TABLE_2.csv")
write_csv_utf8(table2, table2_file)
copy_canonical(table2_file, "HD_T7B_TABLE_2.csv")

s1 <- r4_claims
s1_file <- out_path("06_supplementary_tables", "HD_T7B_SUPPLEMENTARY_TABLE_S1.csv")
write_csv_utf8(s1, s1_file)
copy_canonical(s1_file, "HD_T7B_SUPPLEMENTARY_TABLE_S1.csv")

# -------------------------------------------------------------------------
# Supplementary package assembled from frozen outputs
# -------------------------------------------------------------------------

supp_fig_map <- data.frame(
  supplementary_id = sprintf("Figure S%d", 1:12),
  description = c(
    "Raw library totals", "mRNA PCA", "miRNA PCA from T1", "mRNA quality weights",
    "mRNA adjusted-effect distribution", "mRNA primary versus reduced-zone estimates",
    "mRNA transcript concordance", "Competitive camera versus patient-score module evidence",
    "miRNA quality weights", "miRNA PCA from T3", "miRNA primary versus reduced-zone estimates",
    "miRNA leave-one-patient-out ranges"
  ),
  source_stem = c(
    "results/HD_T1/05_qc_figures/HD_T1_FIG02_RAW_TOTALS",
    "results/HD_T1/05_qc_figures/HD_T1_FIG05_ENST_PCA",
    "results/HD_T1/05_qc_figures/HD_T1_FIG06_MIRNA_PCA",
    "results/HD_T2/08_figures/HD_T2_FIG02_SAMPLE_QUALITY_WEIGHTS",
    "results/HD_T2/08_figures/HD_T2_FIG03_ADJUSTED_EFFECT_DISTRIBUTION",
    "results/HD_T2/08_figures/HD_T2_FIG05_PRIMARY_VS_REDUCED_ZONE",
    "results/HD_T2/08_figures/HD_T2_FIG10_GENE_TRANSCRIPT_CONCORDANCE",
    "results/HD_T2/08_figures/HD_T2_FIG06_PREDEFINED_MODULE_EFFECTS",
    "results/HD_T3/07_figures/HD_T3_FIG03_SAMPLE_QUALITY_WEIGHTS",
    "results/HD_T3/07_figures/HD_T3_FIG04_MIRNA_PCA",
    "results/HD_T3/07_figures/HD_T3_FIG07_PRIMARY_VS_REDUCED_ZONE",
    "results/HD_T3/07_figures/HD_T3_FIG11_LEAVE_ONE_OUT_RANGES"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(supp_fig_map))) {
  for (ext in c(".pdf", ".png")) {
    src <- abs_path(paste0(supp_fig_map$source_stem[i], ext))
    if (!file.exists(src)) stop("Missing supplementary figure source: ", src)
    dest <- out_path("05_supplementary_figures", sprintf("HD_T7B_SUPPLEMENTARY_FIGURE_S%02d%s", i, ext))
    if (!file.copy(src, dest, overwrite = TRUE)) stop("Failed supplementary figure copy: ", src)
  }
}
write_csv_utf8(supp_fig_map, out_path("05_supplementary_figures", "HD_T7B_SUPPLEMENTARY_FIGURE_SOURCE_MAP.csv"))

supp_table_map <- data.frame(
  supplementary_id = c(
    "Table S2", "Table S3a", "Table S3b", "Table S3c", "Table S3d",
    "Table S4a", "Table S4b", "Table S4c", "Table S4d",
    "Table S5a", "Table S5b", "Table S6", "Table S7a", "Table S7b",
    "Table S7c", "Table S7d", "Table S8a", "Table S8b", "Table S9a",
    "Table S9b", "Table S9c", "Table S9d", "Table S10a", "Table S10b"
  ),
  description = c(
    "Patient manifest", "mRNA completely unadjusted", "mRNA no-zone", "mRNA full adjusted", "mRNA reduced-zone",
    "mRNA stability classification", "mRNA transcript concordance", "mRNA warning-sample sensitivity", "mRNA leave-one-out summary",
    "Module camera and patient-score results", "Module stability results", "miRNA filtering and family definitions",
    "Mature-miRNA completely unadjusted", "Mature-miRNA no-zone", "Mature-miRNA full adjusted", "Mature-miRNA reduced-zone",
    "Precursor miRNA adjusted results", "Unresolved miRNA adjusted results", "miRNA stability classification",
    "miRNA zone-dependence audit", "miRNA warning-sample sensitivity", "miRNA leave-one-out summary",
    "mRNA sample influence", "miRNA sample influence"
  ),
  source_file = c(
    "HD_DATASET_AND_SAMPLE_MANIFEST.csv",
    "results/HD_T2/HD_T2_GENE_DE_UNADJUSTED.csv.gz",
    "results/HD_T2/HD_T2_GENE_DE_NO_ZONE_ADJUSTMENT.csv.gz",
    "results/HD_T2/HD_T2_GENE_DE_ADJUSTED.csv.gz",
    "results/HD_T2/HD_T2_GENE_DE_REDUCED_ZONE.csv.gz",
    "results/HD_T2/HD_T2_GENE_STABILITY_CLASSIFICATION.csv",
    "results/HD_T2/HD_T2_GENE_TRANSCRIPT_CONCORDANCE.csv",
    "results/HD_T2/HD_T2_WARNING_SAMPLE_SENSITIVITY.csv",
    "results/HD_T2/HD_T2_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv",
    "results/HD_T2/HD_T2_MODULE_DIFFERENTIAL_RESULTS.csv",
    "results/HD_T2/HD_T2_MODULE_STABILITY_RESULTS.csv",
    "results/HD_T3/02_filtering/HD_T3_FILTER_COUNTS_BY_FAMILY.csv",
    "results/HD_T3/HD_T3_MIRNA_DE_UNADJUSTED.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_DE_NO_ZONE_ADJUSTMENT.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_DE_ADJUSTED.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_DE_REDUCED_ZONE.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_DE_ADJUSTED.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_DE_ADJUSTED.csv.gz",
    "results/HD_T3/HD_T3_MIRNA_STABILITY_CLASSIFICATION.csv",
    "results/HD_T3/HD_T3_ZONE_DEPENDENCE_AUDIT.csv",
    "results/HD_T3/HD_T3_WARNING_SAMPLE_SENSITIVITY.csv",
    "results/HD_T3/HD_T3_LEAVE_ONE_PATIENT_OUT_SUMMARY.csv",
    "results/HD_T2/HD_T2_SAMPLE_INFLUENCE_METRICS.csv",
    "results/HD_T3/HD_T3_SAMPLE_INFLUENCE_METRICS.csv"
  ),
  row_filter = c(
    rep("all rows", 16),
    "feature_family == PRECURSOR",
    "feature_family == UNRESOLVED",
    rep("all rows", 6)
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(supp_table_map))) {
  src <- abs_path(supp_table_map$source_file[i])
  if (!file.exists(src)) stop("Missing supplementary table source: ", src)
  id_safe <- gsub("[^A-Za-z0-9]+", "_", supp_table_map$supplementary_id[i])
  ext <- if (grepl("\\.gz$", src, ignore.case = TRUE)) ".csv.gz" else ".csv"
  dest <- out_path("06_supplementary_tables", paste0("HD_T7B_SUPPLEMENTARY_", id_safe, ext))
  if (supp_table_map$row_filter[i] == "all rows") {
    if (!file.copy(src, dest, overwrite = TRUE)) stop("Failed supplementary table copy: ", src)
  } else {
    d <- read_csv(src)
    fam <- sub("feature_family == ", "", supp_table_map$row_filter[i], fixed = TRUE)
    d <- d[d$feature_family == fam, ]
    dest <- sub("\\.gz$", "", dest)
    write_csv_utf8(d, dest)
  }
}
write_csv_utf8(supp_table_map, out_path("06_supplementary_tables", "HD_T7B_SUPPLEMENTARY_TABLE_SOURCE_MAP.csv"))

# -------------------------------------------------------------------------
# Provenance ledger
# -------------------------------------------------------------------------

hash_sources <- function(keys) {
  inv <- source_inventory[source_inventory$source_key %in% keys, ]
  paste(inv$source_file_sha256[match(keys, inv$source_key)], collapse = ";")
}
source_names <- function(keys) {
  inv <- source_inventory[source_inventory$source_key %in% keys, ]
  paste(inv$source_file[match(keys, inv$source_key)], collapse = ";")
}

provenance <- data.frame(
  figure_id = c(rep("Figure 1", 4), rep("Figure 2", 4), rep("Figure 3", 4), rep("Figure 4", 4), "Table 1", "Table 2", "Supplementary Table S1"),
  panel_id = c(LETTERS[1:4], LETTERS[1:4], LETTERS[1:4], LETTERS[1:4], "NA", "NA", "NA"),
  scientific_question = c(
    "What was measured and paired?", "How was the zone class derived?", "How are disease groups distributed across zones?", "Which model contexts are distinct?",
    "How many genes are FDR significant in each model?", "How do no-zone and full-adjusted gene effects compare?", "How often do no-zone gene signals attenuate or reverse?", "What stability and patient-influence patterns limit gene claims?",
    "How were miRNA families filtered and how do mature-miRNA FDR counts migrate?", "How do no-zone and full-adjusted mature-miRNA effects compare?", "How often do no-zone mature-miRNA signals attenuate or reverse?", "What mature-miRNA stability and patient-influence patterns limit claims?",
    "Do mRNA and mature miRNA show parallel FDR migration?", "Do attenuation and reversal proportions align across analytes?", "How do frozen evidence classes compare without equating their rules?", "What evidence boundary is supported?",
    "What are cohort, covariate and model characteristics?", "What is the cross-model analyte summary?", "Which claims and wording are supported?"
  ),
  source_keys = c(
    "manifest", "r3_coverage", "manifest;t1_zone", "t1_design;t2_reduced;t3_reduced",
    "t2_unadjusted;t2_no_zone;t2_adjusted;t2_reduced", "t2_no_zone;t2_adjusted", "t2_no_zone;t2_adjusted", "t2_stability;t2_loo;t2_influence",
    "t3_filter;t3_unadjusted;t3_no_zone;t3_adjusted;t3_reduced", "t3_no_zone;t3_adjusted", "t3_no_zone;t3_adjusted", "t3_stability;t3_loo;t3_influence",
    "t2_unadjusted;t2_no_zone;t2_adjusted;t2_reduced;t3_unadjusted;t3_no_zone;t3_adjusted;t3_reduced",
    "t2_no_zone;t2_adjusted;t3_no_zone;t3_adjusted",
    "t2_stability;t3_stability",
    "r4_claims;t2_ledger;t3_ledger",
    "manifest;t1_balance;t1_zone", "t2_ledger;t2_stability;t3_ledger;t3_stability", "r4_claims"
  ),
  source_columns = c(
    "patient_id,disease_group,mRNA_present,miRNA_present", "analysis_item,evidence,implication_for_R3", "tissue_zone,disease_group,n", "design-matrix columns;model_name;n_patients",
    "BH_FDR", "ensembl_gene_id,log2_fold_change,BH_FDR", "log2_fold_change,BH_FDR", "stability_class,direction_retention_fraction,maximum_absolute_logFC_change,warning_class",
    "feature_family,before_n,primary_after_n,BH_FDR", "original_miRNA_id,feature_family,logFC,BH_FDR", "logFC,BH_FDR", "feature_family,stability_class,direction_retention_fraction,maximum_absolute_logFC_change,high_disease_effect_influence",
    "BH_FDR,feature_family", "log2_fold_change/logFC,BH_FDR", "stability_class,feature_family", "claim,evidence,allowed_wording,forbidden_wording",
    "sex,BMI,tissue_zone,disease_group", "ledger values;stability_class", "all columns"
  ),
  row_filter = c(
    "all 38 patients", "zone-provenance rows", "all disease-by-zone cells", "frozen model definitions",
    "all tested genes; BH_FDR < 0.05 counts", "all 12,161 genes", "5,389 no-zone FDR-significant genes", "all 12,161 genes and all 38 patient-influence rows",
    "all families for filtering; mature family for FDR ladder", "670 mature miRNAs", "45 no-zone FDR-significant mature miRNAs", "670 mature miRNAs and all 38 patient-influence rows",
    "mRNA genes and mature miRNAs", "no-zone FDR-significant families", "mutually exclusive frozen classes; mature family only for miRNA", "all confirmed claim rows",
    "all 38 patients", "mRNA genes and mature miRNAs", "all 16 claims"
  ),
  model_name = c(
    rep("DESCRIPTIVE", 3), "FOUR_MODEL_LADDER",
    "FOUR_MODEL_LADDER", "NO_ZONE_VS_PRIMARY_ADJUSTED", "NO_ZONE_VS_PRIMARY_ADJUSTED", "PRIMARY_AND_SENSITIVITY",
    "FILTERING_AND_FOUR_MODEL_LADDER", "NO_ZONE_VS_PRIMARY_ADJUSTED", "NO_ZONE_VS_PRIMARY_ADJUSTED", "PRIMARY_AND_SENSITIVITY",
    "FOUR_MODEL_LADDER", "NO_ZONE_VS_PRIMARY_ADJUSTED", "FROZEN_ANALYTE_SPECIFIC_CLASSIFICATION", "EVIDENCE_SYNTHESIS",
    "DESCRIPTIVE", "FOUR_MODEL_LADDER_AND_FROZEN_CLASSIFICATION", "CLAIM_EVIDENCE_MATRIX"
  ),
  statistical_unit = c(
    rep("patient/library", 4), rep("patient/library; one estimate per tested gene", 4),
    rep("patient/library; one estimate per mature miRNA", 4),
    rep("analyte-level summary of patient-level models", 4),
    "patient/library", "patient/library; analyte-level summary", "claim/evidence row"
  ),
  n_total = c(rep(38, 18), 16),
  n_case = c(rep(20, 18), NA),
  n_control = c(rep(18, 18), NA),
  displayed_number = c(
    "38;20;18;12,161;697;670", "three expression-derived zone classes", "2/10;4/5;14/3", "38;38;38;26",
    "6,111;5,389;0;0", "12,161 estimates", "5,255 attenuated;1,016 reversals;0.513→0.149", "4,210 low precision;0 robust;38 LOO summaries",
    "1,806→697;670 mature;46→45→0→0", "670 estimates", "41 attenuated;5 reversals;1.767→0.537", "229 low precision;333 sample-sensitive;93 zone-dependent;0 adjusted-stable",
    "6,111/46;5,389/45;0/0;0/0", "97.5%/91.1% attenuation;18.9%/11.1% reversals", "analyte-specific class percentages", "0 adjusted FDR features;bounded conclusion",
    "38;20/18;sex;BMI;zone;paired analytes", "all cross-model headline counts", "16 claim/evidence rows"
  ),
  transformation_for_plotting = c(
    "counts and diagram layout", "text extraction and diagram layout", "grouped counts", "formula text and model ordering",
    "FDR count by frozen model", "identity join by gene ID", "absolute-effect and sign comparisons already defined by HD-R4", "frozen class counts and ranked patient influence metric",
    "family counts plus mature-family FDR counts", "identity join by miRNA ID", "absolute-effect and sign comparisons already defined by HD-R4", "mature-family frozen class counts and ranked patient influence metric",
    "side-by-side faceted counts", "within-family percentages", "within-family proportions", "text-only evidence ladder",
    "descriptive summaries", "cross-analyte row assembly", "column-preserving copy"
  ),
  interpretation_level = c(
    rep("DESCRIPTIVE_DESIGN_EVIDENCE", 4), rep("FORMAL_FEATURE_MODEL_AND_SENSITIVITY", 4),
    rep("FORMAL_MATURE_FAMILY_AND_SENSITIVITY", 4), rep("INTEGRATED_DESCRIPTIVE_SYNTHESIS", 4),
    "DESCRIPTIVE", "FORMAL_SUMMARY", "CLAIM_BOUNDARY"
  ),
  claim_boundary = c(
    "No directly measured anatomy or causal zone effect", "Expression-derived, not independently pathology confirmed", "Imbalance does not prove causation", "No model is declared uniquely true",
    "6,111 is context; 5,389→0 isolates zone addition", "No disease-gene labels", "Attenuation does not prove zone explains all differences", "Low precision is not robustness or biomarker evidence",
    "697 is all retained records; 670 is the mature formal family", "No independent replication or target inference", "No biomarker miRNA claims", "Mature-only main classification; all-family counts are supplementary",
    "Parallel sensitivity is not regulatory integration", "Within-analyte denominators remain separate", "Classification rules remain analyte specific", "No mechanism, diagnosis or causal adjustment",
    "No healthy-control wording", "No disease-specific signature", "Allowed and prohibited wording preserved"
  ),
  notes = "All displayed values trace to frozen HD-T1–HD-T3 or HD-R4 files; no new model was fitted.",
  stringsAsFactors = FALSE
)

split_keys <- strsplit(provenance$source_keys, ";", fixed = TRUE)
provenance$source_file <- vapply(split_keys, source_names, character(1))
provenance$source_file_sha256 <- vapply(split_keys, hash_sources, character(1))
provenance$source_keys <- NULL
provenance <- provenance[, c(
  "figure_id", "panel_id", "scientific_question", "source_file", "source_file_sha256",
  "source_columns", "row_filter", "model_name", "statistical_unit", "n_total", "n_case",
  "n_control", "displayed_number", "transformation_for_plotting", "interpretation_level",
  "claim_boundary", "notes"
)]
prov_file <- out_path("01_source_provenance", "HD_T7B_PANEL_SOURCE_PROVENANCE.csv")
write_csv_utf8(provenance, prov_file)
copy_canonical(prov_file, "HD_T7B_PANEL_SOURCE_PROVENANCE.csv")

# -------------------------------------------------------------------------
# Numeric audit
# -------------------------------------------------------------------------

audit <- data.frame(
  audit_item = c(
    "patients_total", "patients_case", "patients_control", "tested_mRNA_genes",
    "miRNA_records_before_filter", "miRNA_records_after_filter_all_families",
    "mature_miRNA_after_filter", "mRNA_unadjusted_FDR", "mRNA_no_zone_FDR",
    "mRNA_full_adjusted_FDR", "mRNA_reduced_zone_FDR", "mRNA_attenuated_no_zone_signals",
    "mRNA_direction_reversals", "mRNA_low_precision", "mRNA_robust",
    "mature_miRNA_unadjusted_FDR", "mature_miRNA_no_zone_FDR", "mature_miRNA_full_adjusted_FDR",
    "mature_miRNA_reduced_zone_FDR", "mature_miRNA_attenuated_no_zone_signals",
    "mature_miRNA_direction_reversals", "mature_miRNA_low_precision",
    "mature_miRNA_sample_sensitive_class", "mature_miRNA_zone_dependent_class",
    "mature_miRNA_model_dependent_class", "mature_miRNA_adjusted_stable",
    "all_filtered_miRNA_sample_sensitive_class", "all_filtered_miRNA_zone_dependent_class",
    "stable_camera_modules", "patient_level_module_score_FDR_significant"
  ),
  expected_value = c(
    38, 20, 18, 12161, 1806, 697, 670, 6111, 5389, 0, 0, 5255, 1016, 4210, 0,
    46, 45, 0, 0, 41, 5, 229, 333, 93, 10, 0, 348, 95, 3, 0
  ),
  observed_value = c(
    n_total, n_case, n_control, nrow(t2_adj), sum(t3_filter$before_n), nrow(t3_adj),
    mirna_mature_count(t3_adj), gene_counts[["unadjusted"]], gene_counts[["no_zone"]],
    gene_counts[["full"]], gene_counts[["reduced"]], gene_attenuated, gene_reversed,
    sum(t2_stab$stability_class == "DIRECTIONALLY_STABLE_LOW_PRECISION"),
    sum(t2_stab$stability_class == "ROBUST"), mirna_counts[["unadjusted"]], mirna_counts[["no_zone"]],
    mirna_counts[["full"]], mirna_counts[["reduced"]], mirna_attenuated, mirna_reversed,
    sum(t3_stab_mature$stability_class == "DIRECTIONALLY_STABLE_LOW_PRECISION"),
    sum(t3_stab_mature$stability_class == "SAMPLE_SENSITIVE"),
    sum(t3_stab_mature$stability_class == "ZONE_DEPENDENT"),
    sum(t3_stab_mature$stability_class == "MODEL_DEPENDENT"),
    sum(t3_stab_mature$adjusted_stable),
    sum(t3_stab_filtered$stability_class == "SAMPLE_SENSITIVE"),
    sum(t3_stab_filtered$stability_class == "ZONE_DEPENDENT"),
    stable_camera, score_fdr
  ),
  source_file = c(
    rep("HD_DATASET_AND_SAMPLE_MANIFEST.csv", 3),
    rep("results/HD_T2/HD_T2_GENE_DE_ADJUSTED.csv.gz", 1),
    rep("results/HD_T3/02_filtering/HD_T3_FILTER_COUNTS_BY_FAMILY.csv", 3),
    rep("results/HD_T2 model tables", 4),
    rep("results/HD_T2 no-zone and adjusted model tables", 2),
    rep("results/HD_T2/HD_T2_GENE_STABILITY_CLASSIFICATION.csv", 2),
    rep("results/HD_T3 model tables", 4),
    rep("results/HD_T3 no-zone and adjusted model tables", 2),
    rep("results/HD_T3/HD_T3_MIRNA_STABILITY_CLASSIFICATION.csv", 7),
    rep("results/HD_T2/HD_T2_MODULE_STABILITY_RESULTS.csv", 2)
  ),
  scope_note = c(
    rep("Frozen patient-level cohort", 3),
    "Formal mRNA family", "All deposited miRNA records", "All filtered miRNA families",
    "Formal mature-miRNA family", rep("Formal mRNA model ladder", 4),
    "Among 5,389 no-zone FDR signals", "Among 5,389 no-zone FDR signals",
    "Frozen mutually exclusive gene classes", "Frozen robust flag",
    rep("Formal mature-miRNA model ladder", 4),
    "Among 45 no-zone FDR signals", "Among 45 no-zone FDR signals",
    rep("Mature-miRNA-only main classification", 5),
    rep("All 697 filtered records; supplementary context only", 2),
    "Competitive camera enrichment", "Patient-level module-score model"
  ),
  stringsAsFactors = FALSE
)
audit$status <- ifelse(as.character(audit$expected_value) == as.character(audit$observed_value), "PASS", "FAIL")
if (any(audit$status == "FAIL")) stop("Numeric audit contains failures")
audit_file <- out_path("01_source_provenance", "HD_T7B_NUMERIC_CONSISTENCY_AUDIT.csv")
write_csv_utf8(audit, audit_file)
copy_canonical(audit_file, "HD_T7B_NUMERIC_CONSISTENCY_AUDIT.csv")

# -------------------------------------------------------------------------
# Legends, notes, style, audits and environment
# -------------------------------------------------------------------------

legend1 <- c(
  "# Fig. 1 | Cohort, expression-derived zone provenance and model framework",
  "",
  "**a,** Patient-level design for 38 independent Combo-Seq libraries: 20 hemorrhoidal-disease cases and 18 anal-fissure surgical controls. Each patient contributes paired mRNA and miRNA measurements; 12,161 genes and 670 mature miRNAs enter the formal feature families. **b,** Provenance of the deposited zone class. The original study used keratinocyte/sebocyte expression signatures, ssGSEA and clustering to infer three histological-zone classes; the variable is therefore expression-derived and is not an independently measured anatomical or pathology-confirmed location. **c,** Disease-by-zone patient counts (anoderm 2/10, transition 4/5, intestinal mucosa 14/3 for case/control), showing sparse and asymmetric overlap. **d,** Four distinct frozen model contexts: completely unadjusted; sex/BMI-adjusted without zone; full adjusted including the expression-derived zone; and reduced-zone sensitivity excluding anoderm (n = 26). Statistical unit: patient/library. No inferential P values are used in this descriptive figure. The figure cannot establish anatomical causation or identify a uniquely correct model. Source data are provided in the HD-T7B figure-data files."
)
legend2 <- c(
  "# Fig. 2 | Migration and stability of mRNA disease-effect estimates",
  "",
  "**a,** Counts of genes at Benjamini–Hochberg FDR < 0.05 across the frozen model ladder. The isolated zone-addition comparison is 5,389 sex/BMI-adjusted no-zone signals to 0 full-adjusted signals; 6,111 is completely unadjusted context. **b,** Disease coefficients for all 12,161 tested genes in the no-zone and full-adjusted models. The dashed line is identity; red points reverse direction after adjustment. **c,** Among the 5,389 no-zone FDR signals, 5,255 (97.5%) attenuate in absolute magnitude and 1,016 (18.9%) reverse direction; median absolute log2 fold change changes from 0.513 to 0.149. **d,** Frozen mutually exclusive stability classes and patient leave-one-out influence summaries. Models used limma–voom with quality weights at the patient/library level; FDR was controlled within the gene family. The deposited zone class was inferred from mRNA expression, so attenuation does not prove that anatomy explains all apparent disease differences. No gene remained FDR significant in the full or reduced-zone analyses, and the 4,210 directionally stable signals are low precision rather than robust biomarkers. Source data are provided in the HD-T7B figure-data files."
)
legend3 <- c(
  "# Fig. 3 | Migration and stability of mature-miRNA disease-effect estimates",
  "",
  "**a,** Filtering and family separation of 1,806 deposited miRNA records into 697 retained records (670 mature, 14 precursor and 13 unresolved), followed by the mature-family FDR ladder: 46 completely unadjusted, 45 sex/BMI-adjusted no-zone, 0 full adjusted and 0 reduced-zone signals. **b,** Disease coefficients for all 670 mature miRNAs in the no-zone and full-adjusted models; mature, precursor and unresolved families remain separate, and 5p/3p arms are not merged. **c,** Among the 45 no-zone mature-miRNA FDR signals, 41 (91.1%) attenuate and 5 (11.1%) reverse direction; median absolute log2 fold change changes from 1.767 to 0.537. **d,** Mature-miRNA-only frozen stability classes and patient leave-one-out influence. Models used limma–voom quality weights at the patient/library level with Benjamini–Hochberg correction within the mature family. The zone classification was derived from mRNA expression, and paired miRNA results are not an independent replication or evidence of mRNA–miRNA regulation. The 229 directionally stable signals are low precision; adjusted-stable mature miRNAs = 0. Source data are provided in the HD-T7B figure-data files."
)
legend4 <- c(
  "# Fig. 4 | Cross-analyte inference boundary",
  "",
  "**a,** Side-by-side model migration for mRNA genes and mature miRNAs. Separate y-scales are used because the tested family sizes differ; exact counts are printed. **b,** Within-analyte attenuation and direction-reversal proportions among sex/BMI-adjusted no-zone FDR signals (5,389 genes and 45 mature miRNAs). **c,** Composition of frozen mutually exclusive evidence classes. Gene and miRNA rules were prespecified separately and are not assumed to be identical. **d,** Non-causal evidence ladder: apparent contrasts in models without the expression-derived zone are evaluated after sex/BMI and zone adjustment, reduced-zone analysis and patient-influence checks, leading to an inference-sensitivity conclusion rather than a mechanism or biomarker claim. All analyses use one 38-patient cohort; paired analytes are not independent replications. Zero full-adjusted FDR features do not prove biological equivalence or absence of disease effects. Competitive camera enrichment, target prediction, networks and cell-source analyses are not included in this figure. Source data are provided in the HD-T7B figure-data files."
)
table1_note <- c(
  "# Table 1 | Cohort, covariates and model contexts",
  "",
  "Values are patient counts unless stated otherwise. The comparator is anal fissure surgical controls, not healthy tissue. The deposited histological-zone classification was inferred in the original study from mRNA expression signatures, ssGSEA and clustering and is not independently pathology confirmed. BMI was modeled as a continuous centered covariate; descriptive P values were not used for covariate selection. Age and batch were unavailable in the recovered public metadata. The reduced-zone sensitivity excludes anoderm and retains 18 hemorrhoidal-disease and 8 control patients."
)
table2_note <- c(
  "# Table 2 | Cross-model mRNA and mature-miRNA summary",
  "",
  "FDR counts use Benjamini–Hochberg correction within the prespecified gene or mature-miRNA family. Completely unadjusted counts (6,111 genes; 46 mature miRNAs) are descriptive context. The isolated addition of the expression-derived zone compares the sex/BMI-adjusted no-zone model with the full adjusted model (5,389→0 genes; 45→0 mature miRNAs). Stability counts are mutually exclusive within each analyte; the mature-miRNA row excludes precursor and unresolved records. Directionally stable low-precision signals are not robust biomarkers, and parallel analyte sensitivity does not imply regulatory integration."
)

legend_files <- c(
  "HD_T7B_FIGURE_1_LEGEND.md", "HD_T7B_FIGURE_2_LEGEND.md",
  "HD_T7B_FIGURE_3_LEGEND.md", "HD_T7B_FIGURE_4_LEGEND.md",
  "HD_T7B_TABLE_1_NOTE.md", "HD_T7B_TABLE_2_NOTE.md"
)
legend_content <- list(legend1, legend2, legend3, legend4, table1_note, table2_note)
for (i in seq_along(legend_files)) {
  p <- out_path("07_legends", legend_files[i])
  write_text_utf8(legend_content[[i]], p)
  copy_canonical(p, legend_files[i])
}

style_guide <- c(
  "# HD-T7B figure style guide",
  "",
  "- Backend: R only (`ggplot2`, `gridExtra`, `svglite`, Cairo PDF and R PNG).",
  "- Final assembled width: 183 mm; panel labels are lowercase bold Arial at approximately 8 pt.",
  "- Base typography: Arial, approximately 6–8 pt at final size; white background; no decorative boxes around panels.",
  "- Disease groups: hemorrhoidal disease blue; anal-fissure surgical controls orange.",
  "- Model coding is fixed across figures: completely unadjusted grey; sex/BMI-adjusted no-zone blue; full adjusted rose; reduced-zone purple.",
  "- Directional migration: attenuation teal; reversal coral. Red/green is never the only distinction.",
  "- Stability classes use fixed colors and are always accompanied by text labels or printed counts.",
  "- No model color implies that a model is uniquely correct.",
  "- Axes show complete labels and statistical units. Dense scatter points use transparency; no feature labels are selected for biological interest.",
  "- No 3D effects, gradients, significance-star emphasis, regulatory arrows, networks, mechanism diagrams, ROC plots or biomarker highlighting.",
  "- Primary export: editable SVG and vector PDF. Raster export: PNG at 600 dpi. Individual panel bundles accompany each assembled main figure."
)
write_text_utf8(style_guide, out_path("HD_T7B_FIGURE_STYLE_GUIDE.md"))

input_audit <- c(
  "# HD-T7B input and source audit",
  "",
  "## Result",
  "",
  "`PASS_ALL_REQUIRED_FROZEN_INPUTS_READABLE_AND_HASHED`",
  "",
  sprintf("- Frozen source files enumerated and hashed: %d.", nrow(source_inventory)),
  "- Every main-figure panel, Main Table 1, Main Table 2 and Supplementary Table S1 has a provenance row.",
  "- All headline numbers were recalculated only as counts or plotting transformations from frozen result columns; no model was fitted and no FDR was recalculated.",
  "- HD-T1, HD-T2 and HD-T3 source directories were read only.",
  "- The main mature-miRNA denominator is 670. The 697-record value is retained only for all-family filtering context.",
  "- All-family stability counts (348 sample-sensitive; 95 zone-dependent) are separated from mature-only main counts (333; 93).",
  "",
  "See `HD_T7B_SOURCE_FILE_INVENTORY.csv`, `HD_T7B_PANEL_SOURCE_PROVENANCE.csv` and `HD_T7B_NUMERIC_CONSISTENCY_AUDIT.csv`."
)
write_text_utf8(input_audit, out_path("HD_T7B_INPUT_AND_SOURCE_AUDIT.md"))

boundary_audit <- c(
  "# HD-T7B analysis-boundary audit",
  "",
  "## Status",
  "",
  "`PASS_NO_NEW_STATISTICAL_ANALYSIS`",
  "",
  "- No statistical model was fitted.",
  "- No filtering rule or feature universe was changed.",
  "- No FDR value was recalculated; existing BH_FDR columns were counted only.",
  "- No stability class was regenerated; frozen class labels were tabulated.",
  "- No new enrichment analysis was run.",
  "- No miRNA–mRNA target, integration or regulatory analysis was run.",
  "- No network, biomarker, diagnostic, machine-learning or single-cell analysis was run.",
  "- No HD-T1, HD-T2 or HD-T3 source result was modified.",
  "- No Figure 5, mechanism figure or biological regulatory arrow was created.",
  "- Camera competitive enrichment remains supplementary and is explicitly separated from non-significant patient-level module scores.",
  "- Manuscript prose drafting was not started."
)
write_text_utf8(boundary_audit, out_path("HD_T7B_ANALYSIS_BOUNDARY_AUDIT.md"))

exec_summary <- c(
  "# HD-T7B executive summary",
  "",
  "HD-T7B assembled four main figures, two main tables and a supplementary package exclusively from frozen HD-T1–HD-T3 and HD-R4 outputs. The principal visual comparison is the sex/BMI-adjusted no-zone model versus the full model containing the deposited expression-derived histological-zone classification: 5,389→0 FDR-significant genes and 45→0 FDR-significant mature miRNAs. Completely unadjusted counts (6,111 and 46) are shown only as context.",
  "",
  sprintf("Among no-zone signals, %s/%s genes and %s/%s mature miRNAs attenuate after zone addition; %s genes and %s mature miRNAs reverse direction.", comma(gene_attenuated), comma(5389), comma(mirna_attenuated), comma(45), comma(gene_reversed), comma(mirna_reversed)),
  "",
  "The package preserves the expression-derived nature and possible endogeneity of the zone class, the anal-fissure surgical-control comparator, sparse within-zone overlap and patient influence. It does not claim anatomical causation, disease-specific biomarkers, regulatory integration, mechanism or biological equivalence.",
  "",
  "Visual inspection remains a separate gate. The automated assembly status at this stage is `PENDING_MANUAL_VISUAL_QA`."
)
write_text_utf8(exec_summary, out_path("HD_T7B_EXECUTIVE_SUMMARY.md"))

supp_index <- c(
  "# HD-T7B supplementary index",
  "",
  "## Supplementary Table S1",
  "",
  "- Claim/evidence boundary copied from the frozen HD-R4 claim/evidence matrix.",
  "",
  "## Supplementary figures",
  ""
)
for (i in seq_len(nrow(supp_fig_map))) {
  supp_index <- c(supp_index, sprintf("- **%s:** %s. Frozen source: `%s`.", supp_fig_map$supplementary_id[i], supp_fig_map$description[i], supp_fig_map$source_stem[i]))
}
supp_index <- c(supp_index, "", "## Supplementary tables", "")
for (i in seq_len(nrow(supp_table_map))) {
  supp_index <- c(supp_index, sprintf("- **%s:** %s. Source: `%s`; filter: `%s`.", supp_table_map$supplementary_id[i], supp_table_map$description[i], supp_table_map$source_file[i], supp_table_map$row_filter[i]))
}
supp_index <- c(
  supp_index, "",
  "## Boundary",
  "",
  "Supplementary materials preserve complete frozen results and sensitivity evidence. Unadjusted feature lists are not candidate biomarker tables. Camera results are competitive gene-set tests; patient-level module scores did not pass FDR. No target network, mechanism, ROC, machine-learning model or single-cell map is included."
)
write_text_utf8(supp_index, out_path("HD_T7B_SUPPLEMENTARY_INDEX.md"))

env_lines <- c(
  sprintf("R: %s", R.version.string),
  sprintf("Platform: %s", R.version$platform),
  sprintf("Locale: %s", Sys.getlocale()),
  sprintf("ggplot2: %s", as.character(packageVersion("ggplot2"))),
  sprintf("gridExtra: %s", as.character(packageVersion("gridExtra"))),
  sprintf("svglite: %s", as.character(packageVersion("svglite"))),
  sprintf("data.table: %s", as.character(packageVersion("data.table"))),
  sprintf("digest: %s", as.character(packageVersion("digest"))),
  "patchwork: not installed; not required because gridExtra performs assembly",
  "ragg: not installed; not required because the R PNG device exports 600-dpi PNG",
  "Backend exclusivity: R generated all panels, assembled figures, SVG, PDF and PNG files"
)
write_text_utf8(env_lines, out_path("09_environment", "HD_T7B_ENVIRONMENT_LOCK.txt"))

session_file <- out_path("09_environment", "HD_T7B_SESSION_INFO.txt")
capture.output(sessionInfo(), file = session_file)
script_copy <- out_path("09_environment", "hd_t7b_assemble_figures.R")
file.copy(abs_path("scripts", "hd_t7b_assemble_figures.R"), script_copy, overwrite = TRUE)

visual_qa_pending <- c(
  "# HD-T7B visual QA report",
  "",
  "## Automated export checks",
  "",
  "- Four assembled vector PDFs generated.",
  "- Four assembled editable SVGs generated.",
  "- Four assembled PNGs generated by the R graphics device at 600 dpi.",
  "- Sixteen individual panel bundles generated as PDF, SVG and 600-dpi PNG.",
  "- Fixed model colors and terminology applied across all main figures.",
  "- No Figure 5, network, mechanism arrow, ROC, biomarker label or unadjusted volcano plot generated.",
  "",
  "## Manual visual inspection",
  "",
  "`PENDING_MANUAL_VISUAL_INSPECTION`",
  "",
  "The final gate must not be set to PASS until the four assembled PNGs have been inspected for clipping, overlaps, legibility, model ordering, misleading scales and the distinction between 6,111/46 context and 5,389/45 zone-addition comparisons."
)
write_text_utf8(visual_qa_pending, out_path("HD_T7B_VISUAL_QA_REPORT.md"))

pending_gate <- c(
  "# HD-T7B final gate",
  "",
  "## Current status",
  "",
  "`HD_T7B_REQUIRES_VISUAL_REVISION`",
  "",
  "The numeric, provenance and analysis-boundary checks passed, but the assembled figures require manual visual inspection before the final PASS token can be issued."
)
write_text_utf8(pending_gate, out_path("HD_T7B_FINAL_GATE.md"))

cat("HD-T7B assembly completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n")
cat("Pending manual visual QA.\n")
