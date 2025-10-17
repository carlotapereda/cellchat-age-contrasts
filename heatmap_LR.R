cat("[START] panelC heatmap script\n")

# ---- Logging setup ----
outdir <- "panelC_heatmap_LR_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(outdir, "panelC_heatmap_plots.log")
cat("[INFO] Output dir: ", normalizePath(outdir), "\n", sep = "")
cat("[INFO] Log file  : ", log_file, "\n", sep = "")

# capture stdout (cat/print). Strings are fine here.
sink(log_file, split = TRUE, append = FALSE)

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
geno_cols <- c("E33"  = "#129c4e", "E44" = "#fb9e06")

############################################################
## Helpers
############################################################
# Extract LR table, keep pathway
get_lr_tbl_ct <- function(obj) {
  subsetCommunication(obj) %>%
    transmute(
      lr_pair = .data$interaction_name_2,  # ligand–receptor pair
      pathway = .data$pathway_name,        # pathway for filtering
      source  = .data$source,
      target  = .data$target,
      prob    = .data$prob,
      pval    = .data$pval
    ) %>%
    filter(!is.na(lr_pair), lr_pair != "", is.finite(prob), is.finite(pval))
}

# Sum contributions per (LR, celltype)
rank_celltype_global <- function(obj, label) {
  lr <- get_lr_tbl_ct(obj)
  if (nrow(lr) == 0) {
    return(tibble(name = character(0), celltype = character(0),
                  contribution = numeric(0), pvalues = numeric(0),
                  pathway = character(0), group = character(0)))
  }
  long <- bind_rows(
    lr %>% transmute(lr_pair, pathway, celltype = source, prob, pval),
    lr %>% transmute(lr_pair, pathway, celltype = target, prob, pval)
  )
  contrib <- long %>%
    group_by(lr_pair, pathway, celltype) %>%
    summarise(contribution = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = label) %>%
    rename(name = lr_pair)
}

# Generic builder: age comparison (upper vs lower)
.compare_age_pair <- function(lower_obj, upper_obj, lower_lab, upper_lab,
                              sex_label, geno_label,
                              eps = 1e-6) {
  df <- bind_rows(
    rank_celltype_global(lower_obj, lower_lab),
    rank_celltype_global(upper_obj, upper_lab)
  )
  pair <- df %>%
    group_by(name, pathway, celltype) %>%
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
  pair %>% select(name, pathway, celltype, Sex, Genotype, log2FC)
}

# Build matrix + metadata (Male left, Female right)
.build_heatmap_matrix <- function(df_all) {
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
  
  # Drop rows with fewer than 2 finite values (prevents hclust crash)
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
    # 🔑 Force Male first, Female second
    Sex      = factor(sex,  levels = c("Male","Female")),
    Genotype = factor(geno, levels = c("E33","E44")),
    CellType = factor(cellt, levels = CELLTYPE_ORDER),
    stringsAsFactors = FALSE
  )
  
  # Reorder so Male always left, Female right
  ord <- order(col_meta$Sex, col_meta$Genotype, col_meta$CellType)
  mat      <- mat[, ord, drop = FALSE]
  col_meta <- col_meta[ord, , drop = FALSE]
  
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
  list(col_fun = col_fun, vmax = vmax)
}

############################################################
## Comparison 1: 12 vs 06 months (filtered)
############################################################
df_12_06 <- bind_rows(
  .compare_age_pair(cellchat.FE306, cellchat.FE312, "06Mo", "12Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE406, cellchat.FE412, "06Mo", "12Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME306, cellchat.ME312, "06Mo", "12Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME406, cellchat.ME412, "06Mo", "12Mo", "Male",   "E44")
) %>% filter(pathway %in% c("GABA-A","GABA-B"),
             celltype == "Inhibitory.Neurons") 


hm1 <- .build_heatmap_matrix(df_12_06)
mat1 <- hm1$mat; col_meta1 <- hm1$col_meta; col_labels1 <- hm1$col_labels
top_ann1 <- .build_top_anno(col_meta1)
cf1 <- .make_col_fun(mat1); col_fun1 <- cf1$col_fun; vmax1 <- cf1$vmax

############################################################
## Comparison 2: 18 vs 12 months (filtered)
############################################################
df_18_12 <- bind_rows(
  .compare_age_pair(cellchat.FE312, cellchat.FE318, "12Mo", "18Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE412, cellchat.FE418, "12Mo", "18Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME312, cellchat.ME318, "12Mo", "18Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME412, cellchat.ME418, "12Mo", "18Mo", "Male",   "E44")
) %>% filter(pathway %in% c("GABA-A","GABA-B"),
             celltype == "Inhibitory.Neurons") 

hm2 <- .build_heatmap_matrix(df_18_12)
mat2 <- hm2$mat; col_meta2 <- hm2$col_meta; col_labels2 <- hm2$col_labels
top_ann2 <- .build_top_anno(col_meta2)
cf2 <- .make_col_fun(mat2); col_fun2 <- cf2$col_fun; vmax2 <- cf2$vmax

############################################################
## Comparison 3: 18 vs 06 months (filtered)
############################################################
df_18_06 <- bind_rows(
  .compare_age_pair(cellchat.FE306, cellchat.FE318, "06Mo", "18Mo", "Female", "E33"),
  .compare_age_pair(cellchat.FE406, cellchat.FE418, "06Mo", "18Mo", "Female", "E44"),
  .compare_age_pair(cellchat.ME306, cellchat.ME318, "06Mo", "18Mo", "Male",   "E33"),
  .compare_age_pair(cellchat.ME406, cellchat.ME418, "06Mo", "18Mo", "Male",   "E44")
) %>% filter(pathway %in% c("GABA-A","GABA-B"),
             celltype == "Inhibitory.Neurons") 

hm3 <- .build_heatmap_matrix(df_18_06)
mat3 <- hm3$mat; col_meta3 <- hm3$col_meta; col_labels3 <- hm3$col_labels
top_ann3 <- .build_top_anno(col_meta3)
cf3 <- .make_col_fun(mat3); col_fun3 <- cf3$col_fun; vmax3 <- cf3$vmax


# ---- Align and sort rows alphabetically across all three matrices ----
common_rows <- Reduce(intersect, list(rownames(mat1), rownames(mat2), rownames(mat3)))
if (length(common_rows) == 0) stop("No common LR rows across the three comparisons after filtering.")

# sort alphabetically
common_rows <- sort(common_rows, method = "radix")

# reindex each matrix to the same row set and order
mat1 <- mat1[common_rows, , drop = FALSE]
mat2 <- mat2[common_rows, , drop = FALSE]
mat3 <- mat3[common_rows, , drop = FALSE]

# Get vmax from all 3 matrices and set colors
all_vals <- c(as.numeric(mat1), as.numeric(mat2), as.numeric(mat3))
all_abs  <- abs(all_vals[is.finite(all_vals)])
vmax_all <- if (length(all_abs)) as.numeric(stats::quantile(all_abs, 0.98, na.rm = TRUE)) else 1
if (!is.finite(vmax_all) || vmax_all <= 0) vmax_all <- 1
col_fun_all <- circlize::colorRamp2(c(-vmax_all, 0, vmax_all),
                                    c("#2C3E99", "#F7F7F7", "#B30000"))

# extract pathway
extract_pathway <- function(x) {
  m <- regmatches(x, regexpr("^GABA-[AB]", x))
  ifelse(nzchar(m), m, NA_character_)
}
row_info <- list(
  pathway = factor(extract_pathway(rownames(mat1)), levels = c("GABA-A","GABA-B"))
)

############################################################
## Shorten row labels
############################################################
# Remove "GABA-A-" or "GABA-B-" from individual rownames
row_labels1 <- sub("^GABA-[AB]-\\s*", "", rownames(mat1))
row_labels2 <- sub("^GABA-[AB]-\\s*", "", rownames(mat2))
row_labels3 <- sub("^GABA-[AB]-\\s*", "", rownames(mat3))

############################################################
## Heatmap 1: 12 vs 06 months
############################################################
ht1 <- Heatmap(
  mat1,
  name = "log2FC",
  col  = col_fun_all,
  na_col = "white",
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  row_split = row_info$pathway,     # <-- split rows by pathway (GABA-A / GABA-B labels)
  show_row_names = TRUE,
  row_labels = row_labels1,         # <-- shortened labels
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
  column_title = "12 vs 06 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title  = "log2FC",
    at     = c(-vmax_all, 0, vmax_all),
    labels = round(c(-vmax_all, 0, vmax_all), 2)
  )
)

############################################################
## Heatmap 2: 18 vs 12 months
############################################################
ht2 <- Heatmap(
  mat2,
  name = "log2FC",
  col  = col_fun_all,
  na_col = "white",
  cluster_rows = FALSE,
  row_split = row_info$pathway,     # <-- same split, keep overall labels
  show_row_names = FALSE,           # row names only on ht1
  row_labels = row_labels2,
  show_column_names   = TRUE,
  column_labels       = col_labels2,
  column_names_rot    = 45,
  column_names_gp     = gpar(fontsize = 8),
  column_split        = col_meta2$Sex,
  top_annotation      = top_ann2,
  cluster_columns     = FALSE,
  rect_gp = gpar(col = NA),
  border  = FALSE,
  column_title = "18 vs 12 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = FALSE
)

############################################################
## Heatmap 3: 18 vs 06 months
############################################################
ht3 <- Heatmap(
  mat3,
  name = "log2FC",
  col  = col_fun_all,
  na_col = "white",
  cluster_rows = FALSE,
  row_split = row_info$pathway,     # <-- same split, keep overall labels
  show_row_names = FALSE,
  row_labels = row_labels3,
  show_column_names   = TRUE,
  column_labels       = col_labels3,
  column_names_rot    = 45,
  column_names_gp     = gpar(fontsize = 8),
  column_split        = col_meta3$Sex,
  top_annotation      = top_ann3,
  cluster_columns     = FALSE,
  rect_gp = gpar(col = NA),
  border  = FALSE,
  column_title = "18 vs 06 months",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  show_heatmap_legend = FALSE
)

############################################################
## Draw combined heatmap
############################################################
cat("[SAVE] panelC_combined_heatmap (PNG)\n")
png(file.path(outdir, "panelC_combined_heatmap.png"),
    width = 10, height = 7, units = "in", res = 300)
draw(
  ht1 + ht2 + ht3,
  padding = unit(c(6, 16, 20, 6), "mm"),
  show_annotation_legend = TRUE,
  heatmap_legend_side = "right"
)
dev.off(); cat("[DONE] panelC_combined_heatmap.png\n")


cat("[DONE] All heatmaps generated and saved to ", normalizePath(outdir), "\n", sep = "")

############################################################
# LOG CLOSING 
############################################################

# Close stdout sinks (use a bounded loop to avoid hangs)
for (i in seq_len(10)) {
  if (sink.number() <= 0) break
  try(sink(), silent = TRUE)
}

cat("[DONE] Logs written to: ", log_file, "\n", sep = "")


