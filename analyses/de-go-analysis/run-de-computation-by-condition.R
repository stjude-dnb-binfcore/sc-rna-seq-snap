#!/usr/bin/env Rscript
#################################################################################
# Within-cell-type DE for a single condition metadata column (condition_value1/2/3)
#################################################################################
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(Seurat)
  library(future)
})

#################################################################################
# Set up directories and paths to root_dir and analysis_dir
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript run-de-computation-by-condition.R <condition_value1|condition_value2|condition_value3> [config_path]")
}
output_suffix <- args[[1L]]
config_path <- if (length(args) >= 2L) args[[2L]] else "../../project_parameters.Config.yaml"

suffix_to_col <- function(yaml, suffix) {
  switch(
    suffix,
    condition_value1 = yaml$condition_value1,
    condition_value2 = yaml$condition_value2,
    condition_value3 = yaml$condition_value3,
    stop("Unknown suffix: ", suffix)
  )
}

config_util <- file.path(getwd(), "util", "function-de-go-config.R")
if (!file.exists(config_util)) {
  config_util <- "util/function-de-go-config.R"
}
source(config_util)
yaml <- read_dego_config_yaml(config_path)
source(file.path(
  yaml$root_dir, "analyses", "de-go-analysis", "util", "function-de-global-interaction.R"
))
mode <- normalize_dego_comparison_mode(yaml$dego_condition_comparison_mode)
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis")
if (!mode %in% c("within_celltype", "both")) {
  message("Skipping ", output_suffix, ": mode=", mode, " does not use within-cell-type precompute")
  quit(save = "no", status = 0)
}

condition_col <- suffix_to_col(yaml, output_suffix)
if (is.null(condition_col) || condition_col %in% c("NULL", ".", "")) {
  message("Skipping ", output_suffix, ": column not configured")
  quit(save = "no", status = 0)
}

source(file.path(analysis_dir, "util", "function-de-within-celltype-comparison.R"))

future_globals_value <- as.numeric(yaml$future_globals_value_dego) * 1024^3
options(future.globals.maxSize = future_globals_value, future.rng.onMisuse = "ignore")

parallel_workers <- yaml$dego_parallel_workers
if (is.null(parallel_workers) || is.na(parallel_workers)) {
  parallel_workers <- max(1L, parallelly::availableCores() - 1L)
} else {
  parallel_workers <- as.integer(parallel_workers)
}

min_cells <- as.integer(yaml$min_cells_per_group_dego_module)
if (is.null(min_cells) || is.na(min_cells)) {
  min_cells <- 25L
}

ident1 <- yaml$dego_ident1
ident2 <- yaml$dego_ident2
if (!is.null(ident1) && ident1 %in% c("NULL", ".", "")) ident1 <- NULL
if (!is.null(ident2) && ident2 %in% c("NULL", ".", "")) ident2 <- NULL

data_file <- file.path(root_dir, "analyses", "cell-types-annotation", "results", "cell_types_annotations_all", "seurat_obj_cell_types_annotations_all.rds")

results_dir <- dego_results_de_dir(analysis_dir, mode)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

message("Loading Seurat object for ", output_suffix)
seurat_obj <- readRDS(data_file)
DefaultAssay(seurat_obj) <- yaml$assay_annotation_module

run_within_celltype_de_for_column(
  seurat_obj = seurat_obj,
  condition_col = condition_col,
  output_suffix = output_suffix,
  results_dir = results_dir,
  resolution_list = yaml$resolution_list_find_markers,
  cell_type_col = yaml$cell_type_label_go_module,
  assay = yaml$assay_annotation_module,
  ident1 = ident1,
  ident2 = ident2,
  min_cells_per_group = min_cells,
  parallel_workers = parallel_workers
)
################################################################################################################
