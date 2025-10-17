
cat("[START] panelB heatmap script\n")

# ---- Logging setup ----
outdir <- "panelB_heatmap_pathway_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(outdir, "panelB_heatmap_plots.log")
cat("[INFO] Output dir: ", normalizePath(outdir), "\n", sep = "")
cat("[INFO] Log file  : ", log_file, "\n", sep = "")

# 1) capture stdout (cat/print). Strings are fine here.
sink(log_file, split = TRUE, append = FALSE)

# 2) capture messages/warnings (message(), warnings, package startup msgs).
#    Must use an OPEN CONNECTION for type="message".
msg_con <- file(log_file, open = "at")   # "append text"
sink(msg_con, type = "message")

# Ensure sinks close even if the script errors
on.exit({
  try(sink(type = "message"))  # turn off message sink
  try(close(msg_con))          # close the connection
  try(sink())                  # turn off stdout sink
}, add = TRUE)


# Load CellChat
library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(colorRamp2)
library(CellChat)
cat("[INFO] Libraries loaded (tidyverse, ComplexHeatmap, dplyr, colorRamp2, CellChat)\n")

# Use current working directory (final_code) and only keep files that exist here
order_levels <- c(
  "CC_DBv2_Female_06Mo_E33.rds","CC_DBv2_Female_06Mo_E44.rds",
  "CC_DBv2_Female_12Mo_E33.rds","CC_DBv2_Female_12Mo_E44.rds",
  "CC_DBv2_Female_18Mo_E33.rds","CC_DBv2_Female_18Mo_E44.rds",
  "CC_DBv2_Male_06Mo_E33.rds","CC_DBv2_Male_06Mo_E44.rds",
  "CC_DBv2_Male_12Mo_E33.rds","CC_DBv2_Male_12Mo_E44.rds",
  "CC_DBv2_Male_18Mo_E33.rds","CC_DBv2_Male_18Mo_E44.rds"
)

here_files <- list.files(".", pattern = "^CC_DBv2_.*\\.rds$", full.names = FALSE)
cat("[INFO] Found ", length(here_files), " CC files in cwd: ", paste(here_files, collapse=", "), "\n", sep = "")
if (!length(here_files)) stop("No CC_DBv2_*.rds files found in current dir: ", getwd())

# keep only those in our expected order_levels, preserve that order
keep <- order_levels[order_levels %in% here_files]
cat("[INFO] Using (ordered): ", paste(keep, collapse=", "), "\n", sep = "")
if (!length(keep)) stop("None of the expected files are present in: ", getwd())

files <- keep  # relative names in current directory
cat("[INFO] Current dir: ", getwd(), "\n", sep = "")
cat("[INFO] Will load: ", paste(files, collapse = ", "), "\n", sep = "")

# Function to build object names like cellchat.F06E33, cellchat.M18E44, etc.
make_name <- function(path) {
  # Extract parts from file name
  fname <- basename(path) # e.g. CC_DBv2_Female_06Mo_E33.rds
  fname <- sub("\\.rds$", "", fname)  # remove the .rds extension
  parts <- strsplit(fname, "_")[[1]]
  sex   <- substr(parts[3], 1, 1)   # "F" or "M"
  age   <- gsub("Mo", "", parts[4]) # "06", "12", "18"
  geno  <- parts[5]                 # "E33" or "E44"
  geno_short <- sub("^E(\\d)\\d$", "E\\1", geno)
  paste0("cellchat.", sex, geno_short, age) # NOTE THIS IS DIFFERENT THAT NETDIFFINTERACTION
}

# Load, update, and assign
for (f in files) {
  obj_name <- make_name(f)
  cc <- readRDS(f)
  cc <- updateCellChat(cc)
  assign(obj_name, cc, envir = .GlobalEnv)
}
cat("[INFO] All CellChat objects loaded and updated\n")
cat("[DEBUG] Loaded objects:\n",
    paste(ls(pattern="^cellchat\\."), collapse=", "), "\n", sep = "")


############################################################
## Config
############################################################
cat("[CONFIG] Setting column orders, labels, and colors\n")

# Column panel order for heatmaps
COL_ORDER <- c("Female E33","Female E44","Male E33","Male E44")

# Cell-type ordering + short labels
CELLTYPE_ORDER <- c("Astrocytes","Excitatory.Neurons","Inhibitory.Neurons",
                    "Microglia","Oligodendrocyte.Precursor","Oligodendrocytes")
CELLTYPE_MAP <- c(
  "Astrocytes"                = "Ast",
  "Excitatory.Neurons"        = "Ex.Neu",
  "Inhibitory.Neurons"        = "In.Neu",
  "Microglia"                 = "Mic",
  "Oligodendrocyte.Precursor" = "OPCs",
  "Oligodendrocytes"          = "Oli"
)

# Annotation colors
sex_cols  <- c("Female" = "#D8AADC", "Male" = "#B6D7A8")  # purple/green
geno_cols <- c("E33"  = "#129c4e", "E44" = "#fb9e06")     # blue/red

############################################################
## Helpers
############################################################
cat("[HELPERS] Registering helper functions\n")

# Extract LR table
get_lr_tbl_ct <- function(obj) {
  subsetCommunication(obj) %>%
    transmute(
      pathway = .data$pathway_name,
      source  = .data$source,
      target  = .data$target,
      prob    = .data$prob,
      pval    = .data$pval
    ) %>%
    filter(!is.na(pathway), pathway != "", is.finite(prob), is.finite(pval))
}

# Sum contributions per (pathway, celltype)
rank_celltype_global <- function(obj, label) {
  lr <- get_lr_tbl_ct(obj)
  if (nrow(lr) == 0) {
    cat("[WARN] Empty LR table for object; returning empty tibble\n")
    return(tibble(name = character(0), celltype = character(0),
                  contribution = numeric(0), pvalues = numeric(0), group = character(0)))
  }
  long <- bind_rows(
    lr %>% transmute(pathway, celltype = source, prob, pval, group = label),
    lr %>% transmute(pathway, celltype = target, prob, pval, group = label)
  )
  contrib <- long %>%
    group_by(pathway, celltype, group) %>%
    summarise(contribution = sum(prob, na.rm = TRUE), .groups = "drop")
}

# Generic builder: age comparison (upper vs lower)
.compare_age_pair <- function(lower_obj, upper_obj, lower_lab, upper_lab,
                              sex_label, geno_label,
                              abs_logfc_thr = NA,
                              eps = 1e-6, keep_nonsig_as_na = FALSE) {
  
  df <- bind_rows(
    rank_celltype_global(lower_obj, lower_lab),
    rank_celltype_global(upper_obj, upper_lab)
  ) %>%
    dplyr::rename(name = pathway)  # downstream expects 'name'
  
  cat("[PAIR] ", sex_label, " ", geno_label, " ", upper_lab, " vs ", lower_lab,
      " -> nrows=", nrow(df), "\n", sep = "")
  
  pair <- df %>%
    group_by(name, celltype) %>%
    summarise(
      contrib_low  = sum(contribution[group == lower_lab], na.rm = TRUE),
      contrib_high = sum(contribution[group == upper_lab], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log2FC   = log2((contrib_high + eps) / (contrib_low + eps)),
      Sex      = sex_label,
      Genotype = geno_label
    )
  
  # Optional magnitude threshold (only applied if provided)
  if (!is.na(abs_logfc_thr)) {
    pair <- pair %>%
      mutate(log2FC = ifelse(abs(log2FC) >= abs_logfc_thr, log2FC, NA_real_))
  }
  
  # NOTE: no significance filtering here (no sig_any / p_cut)
  out <- pair %>% select(name, celltype, Sex, Genotype, log2FC)
  cat("[PAIR] Completed ", sex_label, " ", geno_label, " ", upper_lab, " vs ", lower_lab,
      " -> kept=", nrow(out), "\n", sep = "")
  out
}

# Build matrix + metadata with Male left, Female right
.build_heatmap_matrix <- function(df_all) {
  cat("[MATRIX] Input rows: ", nrow(df_all), "\n", sep = "")
  df_all <- df_all %>%
    mutate(Panel = paste(Sex, Genotype),
           Panel = factor(Panel, levels = COL_ORDER),
           celltype = factor(celltype, levels = CELLTYPE_ORDER)) %>%
    mutate(col_key = paste(Panel, as.character(celltype), sep = " | "))
  
  # Pivot to wide format
  mat_wide <- df_all %>%
    select(name, col_key, log2FC) %>%
    tidyr::pivot_wider(names_from = col_key, values_from = log2FC) %>%
    arrange(name)
  
  mat <- as.matrix(mat_wide[, -1, drop = FALSE])
  rownames(mat) <- mat_wide$name
  
  # Replace non-finite with NA
  mat[!is.finite(mat)] <- NA
  
  # Drop rows with fewer than 2 finite values (prevents hclust errors)
  before <- nrow(mat)
  mat <- mat[rowSums(!is.na(mat)) > 1, , drop = FALSE]

  if (nrow(mat) == 0) {
    warning("No valid rows left for heatmap matrix after filtering.")
    return(list(mat = mat, col_meta = NULL, col_labels = NULL))
  }
  
  # Column metadata
  parts <- strsplit(colnames(mat), " \\| ")
  panel <- vapply(parts, `[[`, character(1), 1)
  cellt <- vapply(parts, `[[`, character(1), 2)
  sex   <- ifelse(grepl("^Female", panel), "Female", "Male")
  geno  <- ifelse(grepl("E33$", panel), "E33", "E44")
  
  col_meta <- data.frame(
    Panel    = panel,
    Sex      = factor(sex,  levels = c("Male","Female")),   # explicit order
    Genotype = factor(geno, levels = c("E33","E44")),
    CellType = factor(cellt, levels = CELLTYPE_ORDER),
    stringsAsFactors = FALSE
  )
  
  # Reorder so Female always left, Male right
  ord <- order(col_meta$Sex, col_meta$Genotype, col_meta$CellType)
  mat      <- mat[, ord, drop = FALSE]
  col_meta <- col_meta[ord, , drop = FALSE]
  
  cat("[MATRIX] Final matrix dim: ", nrow(mat), " x ", ncol(mat), "\n", sep = "")
  
  # Column labels (short names)
  short_map <- CELLTYPE_MAP
  missing   <- setdiff(as.character(col_meta$CellType), names(short_map))
  if (length(missing)) short_map <- c(short_map, setNames(missing, missing))
  col_labels <- unname(short_map[as.character(col_meta$CellType)])
  
  list(mat = mat, col_meta = col_meta, col_labels = col_labels)
}

# Top annotation (Sex + Genotype)
.build_top_anno <- function(col_meta) {
  HeatmapAnnotation(
    Sex      = col_meta$Sex,
    Genotype = col_meta$Genotype,
    col = list(Sex = sex_cols, Genotype = geno_cols),
    annotation_name_side = "left",
    simple_anno_size = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 0)
  )
}

# Shared color ramp
.make_col_fun <- function(mat) {
  abs_vals <- abs(mat)
  vmax <- as.numeric(stats::quantile(abs_vals[is.finite(abs_vals)], 0.98, na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  col_fun <- circlize::colorRamp2(c(-vmax, 0, vmax), c("#2C3E99", "#F7F7F7", "#B30000"))
  cat("[COLFUN] vmax (98th pct): ", round(vmax, 3), "\n", sep = "")
  list(col_fun = col_fun, vmax = vmax)
}

############################################################
## Heatmap 1: 12 vs 06 months
############################################################
cat("[STEP] Building df_12_06\n")
df_12_06 <- bind_rows(
  .compare_age_pair(cellchat.FE306, cellchat.FE312, "06Mo", "12Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE406, cellchat.FE412, "06Mo", "12Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME306, cellchat.ME312, "06Mo", "12Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME406, cellchat.ME412, "06Mo", "12Mo", "Male",   "E44")
)
cat("[INFO] df_12_06 rows: ", nrow(df_12_06), "\n", sep = "")
hm1 <- .build_heatmap_matrix(df_12_06)
mat1 <- hm1$mat; col_meta1 <- hm1$col_meta; col_labels1 <- hm1$col_labels
top_ann1 <- .build_top_anno(col_meta1)
cf1 <- .make_col_fun(mat1); col_fun1 <- cf1$col_fun; vmax1 <- cf1$vmax

ht1 <- Heatmap(
  mat1,
  name = "log2FC (12/06)",
  col  = col_fun1,
  na_col = "white",
  cluster_rows = TRUE, show_row_dend = FALSE, clustering_distance_rows = function(x) dist(replace(x, is.na(x), 0)), 
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = unit(45, "mm"),
  show_column_names   = TRUE,
  column_labels       = col_labels1,
  column_names_rot    = 45,
  column_names_gp     = gpar(fontsize = 8),
  column_split        = col_meta1$Sex,
  top_annotation      = top_ann1,
  cluster_columns     = FALSE,
  rect_gp = gpar(col = NA),
  border  = FALSE,
  column_title = "Comparison: 12 vs 06 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title  = "log2FC",
    at     = c(-vmax1, 0, vmax1),
    labels = round(c(-vmax1, 0, vmax1), 2)
  )
)
## ---- SAVE Heatmap 1: 12 vs 06 ----
cat("[SAVE] heatmap_log2FC_12_vs_06 (PNG)\n")
png(file.path(outdir, "heatmap_log2FC_12_vs_06.png"), width = 7, height = 9, units = "in", res = 300)
draw(ht1, padding = unit(c(6, 8, 20, 6), "mm"), show_annotation_legend = TRUE)
dev.off(); cat("[DONE] heatmap_log2FC_12_vs_06.png\n")

############################################################
## Heatmap 2: 18 vs 12 months
############################################################
cat("[STEP] Building df_18_12\n")
df_18_12 <- bind_rows(
  .compare_age_pair(cellchat.FE312, cellchat.FE318, "12Mo", "18Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE412, cellchat.FE418, "12Mo", "18Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME312, cellchat.ME318, "12Mo", "18Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME412, cellchat.ME418, "12Mo", "18Mo", "Male",   "E44")
)
cat("[INFO] df_18_12 rows: ", nrow(df_18_12), "\n", sep = "")
hm2 <- .build_heatmap_matrix(df_18_12)
mat2 <- hm2$mat; col_meta2 <- hm2$col_meta; col_labels2 <- hm2$col_labels
top_ann2 <- .build_top_anno(col_meta2)
cf2 <- .make_col_fun(mat2); col_fun2 <- cf2$col_fun; vmax2 <- cf2$vmax

ht2 <- Heatmap(
  mat2,
  name = "log2FC (18/12)",
  col  = col_fun2,
  na_col = "white",
  cluster_rows = TRUE,show_row_dend = FALSE, clustering_distance_rows = function(x) dist(replace(x, is.na(x), 0)),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = unit(45, "mm"),
  show_column_names   = TRUE,
  column_labels       = col_labels2,
  column_names_rot    = 45,
  column_names_gp     = gpar(fontsize = 8),
  column_split        = col_meta2$Sex,
  top_annotation      = top_ann2,
  cluster_columns     = FALSE,
  rect_gp = gpar(col = NA),
  border  = FALSE,
  column_title = "Comparison: 18 vs 12 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title  = "log2FC",
    at     = c(-vmax2, 0, vmax2),
    labels = round(c(-vmax2, 0, vmax2), 2)
  )
)
## ---- SAVE Heatmap 2: 18 vs 12 ----
cat("[SAVE] heatmap_log2FC_18_vs_12 (PNG)\n")
png(file.path(outdir, "heatmap_log2FC_18_vs_12.png"), width = 7, height = 9, units = "in", res = 300)
draw(ht2, padding = unit(c(6, 8, 20, 6), "mm"), show_annotation_legend = TRUE)
dev.off(); cat("[DONE] heatmap_log2FC_18_vs_12.png\n")

############################################################
## Heatmap 3: 18 vs 06 months
############################################################
cat("[STEP] Building df_18_06\n")
df_18_06 <- bind_rows(
  .compare_age_pair(cellchat.FE306, cellchat.FE318, "06Mo", "18Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE406, cellchat.FE418, "06Mo", "18Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME306, cellchat.ME318, "06Mo", "18Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME406, cellchat.ME418, "06Mo", "18Mo", "Male",   "E44")
)
cat("[INFO] df_18_06 rows: ", nrow(df_18_06), "\n", sep = "")
hm3 <- .build_heatmap_matrix(df_18_06)
mat3 <- hm3$mat; col_meta3 <- hm3$col_meta; col_labels3 <- hm3$col_labels
top_ann3 <- .build_top_anno(col_meta3)
cf3 <- .make_col_fun(mat3); col_fun3 <- cf3$col_fun; vmax3 <- cf3$vmax

ht3 <- Heatmap(
  mat3,
  name = "log2FC (18/06)",
  col  = col_fun3,
  na_col = "white",
  cluster_rows = TRUE,show_row_dend = FALSE, clustering_distance_rows = function(x) dist(replace(x, is.na(x), 0)),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = unit(45, "mm"),
  show_column_names   = TRUE,
  column_labels       = col_labels3,
  column_names_rot    = 45,
  column_names_gp     = gpar(fontsize = 8),
  column_split        = col_meta3$Sex,
  top_annotation      = top_ann3,
  cluster_columns     = FALSE,
  rect_gp = gpar(col = NA),
  border  = FALSE,
  column_title = "Comparison: 18 vs 06 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title  = "log2FC",
    at     = c(-vmax3, 0, vmax3),
    labels = round(c(-vmax3, 0, vmax3), 2)
  )
)
## ---- SAVE Heatmap 3: 18 vs 06 ----
cat("[SAVE] heatmap_log2FC_18_vs_06 (PNG)\n")
png(file.path(outdir, "heatmap_log2FC_18_vs_06.png"), width = 7, height = 9, units = "in", res = 300)
draw(ht3, padding = unit(c(6, 8, 20, 6), "mm"), show_annotation_legend = TRUE)
dev.off(); cat("[DONE] heatmap_log2FC_18_vs_06.png\n")

cat("[DONE] All heatmaps generated and saved to ", normalizePath(outdir), "\n", sep = "")
# Close message sink(s) if any
while (sink.number(type = "message") > 0) sink(type = "message")
# Close stdout sink(s) if any
while (sink.number() > 0) sink()

# If you opened an explicit connection for messages, close it if still open
if (exists("msg_con") && isOpen(msg_con)) close(msg_con)

cat("[DONE] Logs written to: ", log_file, "\n", sep = "")

# --- ensure devices and process exit cleanly ---
if (!is.null(dev.list())) while (!is.null(dev.list())) dev.off()
flush.console()
quit(save = "no", status = 0, runLast = FALSE)

