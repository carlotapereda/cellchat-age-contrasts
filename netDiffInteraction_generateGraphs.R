#Generate Panel A graphs using netVisual_DiffInteraction function in CellChat
# Function was edited for visual purposes

## =========================
## LOAD FILES
## =========================

# ---- Logging setup ----
outdir <- "panelA_diff_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(outdir, "panelA_diff_plots.log")

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
library(CellChat)
cat("[INFO] Loaded CellChat\n")

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
if (!length(here_files)) stop("No CC_DBv2_*.rds files found in current dir: ", getwd())

# keep only those in our expected order_levels, preserve that order
keep <- order_levels[order_levels %in% here_files]
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
  paste0("cellchat.", sex, age, geno)
}

# Load, update, and assign
for (f in files) {
  obj_name <- make_name(f)
  cat("[LOAD] ", obj_name, " <- ", basename(f), "\n", sep = "")
  cc <- readRDS(f)
  cc <- updateCellChat(cc)
  assign(obj_name, cc, envir = .GlobalEnv)
}
cat("[INFO] All CellChat objects loaded and updated\n")

source("netDiffinteraction_adapted.R")
cat("[INFO] Sourced netDiffinteraction_adapted.R\n")

## =========================
## Merge by Sex × Genotype
## (06, 12, 18 included)
## =========================

## MALE · E33  (ME3)
object.list.ME3 <- list(
  M06E33 = cellchat.M06E33,
  M12E33 = cellchat.M12E33,
  M18E33 = cellchat.M18E33
)
cat("[MERGE] Building ME3\n")
cellchat.ME3 <- mergeCellChat(object.list.ME3, add.names = names(object.list.ME3))
cat("[MERGE] Done ME3\n")

## FEMALE · E33 (FE3)
object.list.FE3 <- list(
  F06E33 = cellchat.F06E33,
  F12E33 = cellchat.F12E33,
  F18E33 = cellchat.F18E33
)
cat("[MERGE] Building FE3\n")
cellchat.FE3 <- mergeCellChat(object.list.FE3, add.names = names(object.list.FE3))
cat("[MERGE] Done FE3\n")

## MALE · E44  (ME4)
object.list.ME4 <- list(
  M06E44 = cellchat.M06E44,
  M12E44 = cellchat.M12E44,
  M18E44 = cellchat.M18E44
)
cat("[MERGE] Building ME4\n")
cellchat.ME4 <- mergeCellChat(object.list.ME4, add.names = names(object.list.ME4))
cat("[MERGE] Done ME4\n")

## FEMALE · E44 (FE4)
object.list.FE4 <- list(
  F06E44 = cellchat.F06E44,
  F12E44 = cellchat.F12E44,
  F18E44 = cellchat.F18E44
)
cat("[MERGE] Building FE4\n")
cellchat.FE4 <- mergeCellChat(object.list.FE4, add.names = names(object.list.FE4))
cat("[MERGE] Done FE4\n")

rm(object.list.ME3, object.list.FE3, object.list.ME4, object.list.FE4)
cat("[CLEANUP] Removed temporary object lists\n")

## Helper function to keep parameters consistent
plot_diff <- function(obj, comp) {
  netVisual_diffInteraction_adapted(
    obj,
    comparison = comp,
    weight.scale = TRUE,
    measure = "weight",
    arrow.size = 0.5,
    vertex.label.cex = 0.9,
    alpha.edge = 1,
    top = 0.9
  )
}

dir.create("panelA_diff_plots", showWarnings = FALSE)
cat("[INFO] Output directory ready: panelA_diff_plots\n")

WIDTH = 1200
HEIGHT = 1200

## MALE E33
cat("[PLOT] ME3: 06 vs 12 -> netDiff_ME3_M06E33_vs_M12E33.png\n")
png("panelA_diff_plots/netDiff_ME3_M06E33_vs_M12E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME3, c(1,2))
dev.off(); cat("[SAVED] netDiff_ME3_M06E33_vs_M12E33.png\n")

cat("[PLOT] ME3: 12 vs 18 -> netDiff_ME3_M12E33_vs_M18E33.png\n")
png("panelA_diff_plots/netDiff_ME3_M12E33_vs_M18E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME3, c(2,3))
dev.off(); cat("[SAVED] netDiff_ME3_M12E33_vs_M18E33.png\n")

cat("[PLOT] ME3: 06 vs 18 -> netDiff_ME3_M06E33_vs_M18E33.png\n")
png("panelA_diff_plots/netDiff_ME3_M06E33_vs_M18E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME3, c(1,3))
dev.off(); cat("[SAVED] netDiff_ME3_M06E33_vs_M18E33.png\n")

## FEMALE E33
cat("[PLOT] FE3: 06 vs 12 -> netDiff_FE3_F06E33_vs_F12E33.png\n")
png("panelA_diff_plots/netDiff_FE3_F06E33_vs_F12E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE3, c(1,2))
dev.off(); cat("[SAVED] netDiff_FE3_F06E33_vs_F12E33.png\n")

cat("[PLOT] FE3: 12 vs 18 -> netDiff_FE3_F12E33_vs_F18E33.png\n")
png("panelA_diff_plots/netDiff_FE3_F12E33_vs_F18E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE3, c(2,3))
dev.off(); cat("[SAVED] netDiff_FE3_F12E33_vs_F18E33.png\n")

cat("[PLOT] FE3: 06 vs 18 -> netDiff_FE3_F06E33_vs_F18E33.png\n")
png("panelA_diff_plots/netDiff_FE3_F06E33_vs_F18E33.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE3, c(1,3))
dev.off(); cat("[SAVED] netDiff_FE3_F06E33_vs_F18E33.png\n")

## MALE E44
cat("[PLOT] ME4: 06 vs 12 -> netDiff_ME4_M06E44_vs_M12E44.png\n")
png("panelA_diff_plots/netDiff_ME4_M06E44_vs_M12E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME4, c(1,2))
dev.off(); cat("[SAVED] netDiff_ME4_M06E44_vs_M12E44.png\n")

cat("[PLOT] ME4: 12 vs 18 -> netDiff_ME4_M12E44_vs_M18E44.png\n")
png("panelA_diff_plots/netDiff_ME4_M12E44_vs_M18E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME4, c(2,3))
dev.off(); cat("[SAVED] netDiff_ME4_M12E44_vs_M18E44.png\n")

cat("[PLOT] ME4: 06 vs 18 -> netDiff_ME4_M06E44_vs_M18E44.png\n")
png("panelA_diff_plots/netDiff_ME4_M06E44_vs_M18E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.ME4, c(1,3))
dev.off(); cat("[SAVED] netDiff_ME4_M06E44_vs_M18E44.png\n")

## FEMALE E44
cat("[PLOT] FE4: 06 vs 12 -> netDiff_FE4_F06E44_vs_F12E44.png\n")
png("panelA_diff_plots/netDiff_FE4_F06E44_vs_F12E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE4, c(1,2))
dev.off(); cat("[SAVED] netDiff_FE4_F06E44_vs_F12E44.png\n")

cat("[PLOT] FE4: 12 vs 18 -> netDiff_FE4_F12E44_vs_F18E44.png\n")
png("panelA_diff_plots/netDiff_FE4_F12E44_vs_F18E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE4, c(2,3))
dev.off(); cat("[SAVED] netDiff_FE4_F12E44_vs_F18E44.png\n")

cat("[PLOT] FE4: 06 vs 18 -> netDiff_FE4_F06E44_vs_F18E44.png\n")
png("panelA_diff_plots/netDiff_FE4_F06E44_vs_F18E44.png", width=WIDTH, height=HEIGHT, res=300)
plot_diff(cellchat.FE4, c(1,3))
dev.off(); cat("[SAVED] netDiff_FE4_F06E44_vs_F18E44.png\n")

cat("[DONE] All Panel A plots generated and saved\n")

# ---- Close logs ----
sink(type = "message")  # stop capturing messages
sink()                  # stop capturing output
cat("[DONE] Logs written to: ", log_file, "\n", sep = "")

