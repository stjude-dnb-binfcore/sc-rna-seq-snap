#!/usr/bin/env Rscript
#################################################################################
# Optional DE pre-computation for de-go-analysis (within_celltype / both modes)
#################################################################################
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(Seurat)
  library(future)
})

#################################################################################
# load config file
run_de_go_computation <- function(config_path = "../../project_parameters.Config.yaml",
                                  module_dir = ".") {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(module_dir)

  config_util <- file.path(module_dir, "util", "function-de-go-config.R")
  if (!file.exists(config_util)) {
    config_util <- "util/function-de-go-config.R"
  }
  source(config_util)
  yaml <- read_dego_config_yaml(config_path)

  root_dir <- yaml$root_dir
  analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis")
  source(file.path(analysis_dir, "util", "function-de-global-interaction.R"))
  source(file.path(analysis_dir, "util", "function-de-within-celltype-comparison.R"))

  mode <- normalize_dego_comparison_mode(yaml$dego_condition_comparison_mode)
  if (!run_dego_mode_needs_precompute(mode)) {
    message("dego_condition_comparison_mode=", mode, ": skipping run-de-computation.R (DE runs in 01-de-analysis.Rmd)")
    return(invisible(FALSE))
  }

  data_file <- file.path(root_dir, "analyses", "cell-types-annotation", "results", "cell_types_annotations_all", "seurat_obj_cell_types_annotations_all.rds")
  
  results_dir <- dego_results_de_dir(analysis_dir, mode)
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  future_globals_value <- as.numeric(yaml$future_globals_value_dego) * 1024^3
  options(future.globals.maxSize = future_globals_value, future.rng.onMisuse = "ignore")
  future::plan(future::multisession, workers = max(1L, parallelly::availableCores() - 1L))

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

  assay <- yaml$assay_annotation_module
  cell_type_col <- yaml$cell_type_label_go_module
  resolution_list <- yaml$resolution_list_find_markers
  ident1 <- yaml$dego_ident1
  ident2 <- yaml$dego_ident2
  if (!is.null(ident1) && ident1 %in% c("NULL", ".", "")) ident1 <- NULL
  if (!is.null(ident2) && ident2 %in% c("NULL", ".", "")) ident2 <- NULL

  message("Loading Seurat object: ", data_file)
  seurat_obj <- readRDS(data_file)
  DefaultAssay(seurat_obj) <- assay

  message("FindAllMarkers: all cell types")
  Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
  all.markers <- FindAllMarkers(
    seurat_obj,
    assay = assay,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    min.pct = 0.1,
    return.thresh = 0.01
  )
  finalize_overall_celltype_markers(
    markers = all.markers,
    idents_counts = table(Idents(seurat_obj)),
    min_cells_per_group = min_cells,
    results_dir = results_dir,
    resolution_list = resolution_list
  )

  condition_jobs <- list(
    list(col = yaml$condition_value1, suffix = "condition_value1"),
    list(col = yaml$condition_value2, suffix = "condition_value2"),
    list(col = yaml$condition_value3, suffix = "condition_value3")
  )
  condition_jobs <- Filter(function(j) {
    !is.null(j$col) && !j$col %in% c("NULL", ".", "") && nzchar(as.character(j$col))
  }, condition_jobs)

  run_within <- mode %in% c("within_celltype", "both")
  run_global <- identical(mode, "both")

  run_parallel_columns <- identical(toupper(as.character(yaml$dego_parallel_condition_columns)), "YES")
  n_col_jobs <- length(condition_jobs)

  run_within_column <- function(job) {
    run_within_celltype_de_for_column(
      seurat_obj = seurat_obj,
      condition_col = job$col,
      output_suffix = job$suffix,
      results_dir = results_dir,
      resolution_list = resolution_list,
      cell_type_col = cell_type_col,
      assay = assay,
      ident1 = ident1,
      ident2 = ident2,
      min_cells_per_group = min_cells,
      parallel_workers = parallel_workers
    )
  }

  if (run_within && run_parallel_columns && n_col_jobs > 1L) {
    message(
      "Running ", n_col_jobs,
      " within-cell-type condition columns in parallel (separate Rscript processes)"
    )
    config_abs <- normalizePath(config_path, mustWork = TRUE)
    module_abs <- normalizePath(module_dir, mustWork = TRUE)
    rscript_path <- file.path(module_abs, "run-de-computation-by-condition.R")
    bg_cmds <- vapply(condition_jobs, function(job) {
      paste0(
        "Rscript --vanilla ",
        shQuote(rscript_path), " ",
        shQuote(job$suffix), " ",
        shQuote(config_abs),
        " > ",
        shQuote(file.path(results_dir, paste0("de_compute_", job$suffix, ".out"))),
        " 2> ",
        shQuote(file.path(results_dir, paste0("de_compute_", job$suffix, ".err"))),
        " &"
      )
    }, character(1L))
    shell_cmd <- paste(c(paste(bg_cmds, collapse = " "), "wait"), collapse = " ")
    exit_status <- system(shell_cmd)
    if (exit_status != 0L) {
      warning("Parallel within-cell-type DE jobs exited with status ", exit_status)
    }
  } else if (run_within) {
    for (job in condition_jobs) {
      run_within_column(job)
    }
  }

  if (run_global) {
    message("FindAllMarkers: global cell_type x condition interaction (both mode)")
    for (job in condition_jobs) {
      run_global_interaction_de_for_column(
        seurat_obj = seurat_obj,
        cell_type_col = cell_type_col,
        condition_col = job$col,
        output_suffix = job$suffix,
        results_dir = results_dir,
        resolution_list = resolution_list,
        assay = assay,
        min_cells_per_group = min_cells
      )
    }
  }

  future::plan(future::sequential)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) >= 1L) args[[1L]] else "../../project_parameters.Config.yaml"
  run_de_go_computation(config_path = config_path, module_dir = ".")
}
