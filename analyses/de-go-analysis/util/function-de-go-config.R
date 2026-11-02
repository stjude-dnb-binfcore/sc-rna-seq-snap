############################################################################################################
#' Resolve de-go-analysis parameters when Rmd is knit without run-de-go-analysis.R
############################################################################################################

is_dego_param_placeholder <- function(x) {
  is.null(x) ||
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(trimws(as.character(x))) ||
    identical(as.character(x), ".")
}

find_project_config_path <- function(root_dir = ".") {
  candidates <- unique(c(
    file.path(root_dir, "project_parameters.Config.yaml"),
    "../../project_parameters.Config.yaml",
    "../../../project_parameters.Config.yaml"
  ))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    return(NULL)
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_dego_project_yaml <- function(root_dir = ".") {
  cfg_path <- find_project_config_path(root_dir)
  if (is.null(cfg_path)) {
    return(NULL)
  }
  read_dego_config_yaml(cfg_path)
}

read_dego_config_yaml <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  lines <- readLines(config_path, encoding = "UTF-8", warn = FALSE)
  if (length(lines) == 0L) {
    stop("Configuration file is empty: ", config_path)
  }
  yaml::read_yaml(text = paste(lines, collapse = "\n"))
}

resolve_dego_string_param <- function(param_value,
                                      yaml_key,
                                      root_dir,
                                      label = yaml_key) {
  if (!is_dego_param_placeholder(param_value)) {
    return(as.character(param_value))
  }
  yaml_proj <- read_dego_project_yaml(root_dir)
  if (!is.null(yaml_proj) && !is.null(yaml_proj[[yaml_key]])) {
    val <- yaml_proj[[yaml_key]]
    if (!is_dego_param_placeholder(val)) {
      return(as.character(val))
    }
  }
  stop(
    label, " is not set. Use bash run-de-go-analysis.sh, or set ", yaml_key,
    " in project_parameters.Config.yaml (and pass root_dir in report params)."
  )
}

normalize_dego_mode_label <- function(mode) {
  mode <- tolower(gsub("-", "_", as.character(mode)))
  valid <- c("global_interaction", "within_celltype", "both")
  if (!mode %in% valid) {
    stop(
      "dego_condition_comparison_mode must be one of: ",
      paste(valid, collapse = ", "),
      "; got: ", mode
    )
  }
  mode
}

#' Mode-specific results: results/results_<mode>/01-de-analysis/
dego_results_de_dir <- function(analysis_dir, mode) {
  mode <- normalize_dego_mode_label(mode)
  file.path(analysis_dir, "results", paste0("results_", mode), "01-de-analysis")
}

#' Mode-specific plots root: plots/plots_<mode>/
dego_module_plots_dir <- function(analysis_dir, mode) {
  mode <- normalize_dego_mode_label(mode)
  file.path(analysis_dir, "plots", paste0("plots_", mode))
}

resolve_dego_min_cells_param <- function(param_value, root_dir = ".") {
  if (!is_dego_param_placeholder(param_value)) {
    val <- suppressWarnings(as.integer(param_value))
    if (length(val) == 1L && is.finite(val) && val > 0L) {
      return(val)
    }
  }
  yaml_proj <- read_dego_project_yaml(root_dir)
  if (!is.null(yaml_proj$min_cells_per_group_dego_module)) {
    val <- suppressWarnings(as.integer(yaml_proj$min_cells_per_group_dego_module))
    if (length(val) == 1L && is.finite(val) && val > 0L) {
      return(val)
    }
  }
  25L
}

add_cluster_cell_counts_to_markers <- function(markers, idents_counts) {
  if (is.null(markers) || nrow(markers) == 0L || !"cluster" %in% names(markers)) {
    return(markers)
  }
  markers$cluster_cell_count <- as.integer(idents_counts[as.character(markers$cluster)])
  markers$cluster_cell_count[is.na(markers$cluster_cell_count)] <- 0L
  markers
}

filter_markers_min_cells <- function(markers, min_cells_per_group) {
  min_cells_per_group <- as.integer(min_cells_per_group)
  if (length(min_cells_per_group) != 1L || !is.finite(min_cells_per_group) || min_cells_per_group <= 0L) {
    return(markers)
  }
  if (is.null(markers) || nrow(markers) == 0L) {
    return(markers)
  }
  if (all(c("n_ident1", "n_ident2") %in% names(markers))) {
    keep <- markers$n_ident1 >= min_cells_per_group & markers$n_ident2 >= min_cells_per_group
    return(markers[keep, , drop = FALSE])
  }
  if ("cluster_cell_count" %in% names(markers)) {
    return(markers[markers$cluster_cell_count >= min_cells_per_group, , drop = FALSE])
  }
  markers
}

filter_comparison_markers_min_cells <- function(markers, min_cells_per_group) {
  filter_markers_min_cells(markers, min_cells_per_group)
}

dego_low_cell_exclusion_path <- function(results_dir, resolution_list, label) {
  file.path(results_dir, paste0("Res_", resolution_list, "_Low_cell_", label, ".tsv"))
}

write_low_cell_exclusion_tsv <- function(markers,
                                         min_cells_per_group,
                                         results_dir,
                                         resolution_list,
                                         label) {
  if (!"cluster_cell_count" %in% names(markers)) {
    return(invisible(NULL))
  }
  low_clusters <- unique(as.character(markers$cluster[markers$cluster_cell_count < min_cells_per_group]))
  if (length(low_clusters) == 0L) {
    return(invisible(NULL))
  }
  excl <- unique(markers[markers$cluster %in% low_clusters, c("cluster", "cluster_cell_count"), drop = FALSE])
  excl$min_cells_per_group <- min_cells_per_group
  excl$reason <- "below_min_cells_per_group"
  out_path <- dego_low_cell_exclusion_path(results_dir, resolution_list, label)
  write.table(excl, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(out_path)
}

finalize_markers_with_min_cells <- function(markers,
                                            idents_counts,
                                            min_cells_per_group,
                                            results_dir,
                                            resolution_list,
                                            exclusion_label,
                                            markers_filename) {
  min_cells_per_group <- resolve_dego_min_cells_param(min_cells_per_group)
  markers <- add_cluster_cell_counts_to_markers(markers, idents_counts)
  excluded_path <- write_low_cell_exclusion_tsv(
    markers = markers,
    min_cells_per_group = min_cells_per_group,
    results_dir = results_dir,
    resolution_list = resolution_list,
    label = exclusion_label
  )
  markers <- filter_markers_min_cells(markers, min_cells_per_group)
  out_path <- file.path(results_dir, markers_filename)
  write.table(markers, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(list(
    markers = markers,
    path = out_path,
    excluded_path = excluded_path,
    min_cells_per_group = min_cells_per_group
  ))
}

finalize_overall_celltype_markers <- function(markers,
                                              idents_counts,
                                              min_cells_per_group,
                                              results_dir,
                                              resolution_list) {
  finalize_markers_with_min_cells(
    markers = markers,
    idents_counts = idents_counts,
    min_cells_per_group = min_cells_per_group,
    results_dir = results_dir,
    resolution_list = resolution_list,
    exclusion_label = "overall",
    markers_filename = paste0("Res_", resolution_list, "_Markers_all.tsv")
  )
}

finalize_global_interaction_markers <- function(markers,
                                                idents_counts,
                                                min_cells_per_group,
                                                results_dir,
                                                resolution_list,
                                                output_suffix) {
  finalize_markers_with_min_cells(
    markers = markers,
    idents_counts = idents_counts,
    min_cells_per_group = min_cells_per_group,
    results_dir = results_dir,
    resolution_list = resolution_list,
    exclusion_label = paste0("global_interaction_", output_suffix),
    markers_filename = paste0("Res_", resolution_list, "_Markers_all_", output_suffix, ".tsv")
  )
}

read_dego_markers_filtered <- function(path, min_cells_per_group) {
  if (!file.exists(path)) {
    stop("Marker file not found: ", path)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  filter_markers_min_cells(df, min_cells_per_group)
}

read_low_cell_exclusion_table <- function(results_dir, resolution_list, label) {
  path <- dego_low_cell_exclusion_path(results_dir, resolution_list, label)
  if (!file.exists(path)) {
    return(NULL)
  }
  read.delim(path, stringsAsFactors = FALSE)
}

dego_within_contrast_path <- function(results_dir, resolution_list, suffix) {
  file.path(results_dir, paste0("Res_", resolution_list, "_Within_contrast_", suffix, ".tsv"))
}

read_within_celltype_contrast <- function(results_dir,
                                          resolution_list,
                                          suffix,
                                          markers = NULL) {
  if (!is.null(markers) && nrow(markers) > 0L && all(c("ident1", "ident2") %in% names(markers))) {
    return(data.frame(
      condition_column = if ("condition_column" %in% names(markers)) {
        as.character(markers$condition_column[[1L]])
      } else {
        NA_character_
      },
      ident1 = as.character(markers$ident1[[1L]]),
      ident2 = as.character(markers$ident2[[1L]]),
      ident1_role = "ident.1 (positive avg_log2FC = higher expression vs ident.2)",
      ident2_role = "ident.2 (reference)",
      all_levels_in_column = NA_character_,
      min_cells_per_group = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  path <- dego_within_contrast_path(results_dir, resolution_list, suffix)
  if (!file.exists(path)) {
    return(NULL)
  }
  read.delim(path, stringsAsFactors = FALSE)
}

within_celltype_contrast_cell_counts <- function(markers) {
  if (is.null(markers) || nrow(markers) == 0L || !"cluster" %in% names(markers)) {
    return(NULL)
  }
  cols <- intersect(c("cluster", "ident1", "ident2", "n_ident1", "n_ident2"), names(markers))
  if (!all(c("cluster", "n_ident1", "n_ident2") %in% cols)) {
    return(NULL)
  }
  unique(markers[, cols, drop = FALSE])
}

cat_within_celltype_contrast_report <- function(contrast,
                                                markers = NULL,
                                                condition_label = NULL) {
  if (is.null(contrast) || nrow(contrast) == 0L) {
    cat("\n\n*Within-cell-type contrast metadata not found. Re-run DE precompute to write `Res_*_Within_contrast_*.tsv`.*\n\n")
    return(invisible(NULL))
  }
  cond_col <- contrast$condition_column[[1L]]
  if (!is.na(condition_label) && nzchar(as.character(condition_label))) {
    cat("\n\n**Condition metadata column:** `", condition_label, "` (`", cond_col, "`)\n\n", sep = "")
  } else {
    cat("\n\n**Condition metadata column:** `", cond_col, "`\n\n", sep = "")
  }
  cat(
    "**FindMarkers contrast (same for every cell type in this section):**  \n",
    "- **ident.1 = ", contrast$ident1[[1L]], "** — genes with positive `avg_log2FC` are higher in this group.  \n",
    "- **ident.2 = ", contrast$ident2[[1L]], "** — reference group.  \n\n",
    sep = ""
  )
  if ("all_levels_in_column" %in% names(contrast) &&
      !is.na(contrast$all_levels_in_column[[1L]]) &&
      nzchar(as.character(contrast$all_levels_in_column[[1L]]))) {
    cat("All levels in column: ", contrast$all_levels_in_column[[1L]], ".\n\n", sep = "")
  }
  counts_tbl <- within_celltype_contrast_cell_counts(markers)
  if (!is.null(counts_tbl) && nrow(counts_tbl) > 0L) {
    cat("**Cells per cell type used in the contrast:**\n\n")
    print(knitr::kable(counts_tbl))
    cat("\n")
  }
  invisible(NULL)
}

summarize_dego_min_cells_exclusions <- function(results_dir, resolution_list) {
  low_cell_pat <- paste0("^Res_", resolution_list, "_Low_cell_.*\\.tsv$")
  skipped_pat <- paste0("^Res_", resolution_list, "_Skipped_celltypes_.*\\.tsv$")
  low_files <- list.files(results_dir, pattern = low_cell_pat, full.names = FALSE)
  skipped_files <- list.files(results_dir, pattern = skipped_pat, full.names = FALSE)
  rows <- list()
  for (f in low_files) {
    label <- sub(paste0("^Res_", resolution_list, "_Low_cell_(.*)\\.tsv$"), "\\1", f)
    tbl <- read.delim(file.path(results_dir, f), stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      analysis = label,
      exclusion_type = "low_cell_count",
      n_excluded = nrow(tbl),
      stringsAsFactors = FALSE
    )
  }
  for (f in skipped_files) {
    label <- sub(paste0("^Res_", resolution_list, "_Skipped_(.*)\\.tsv$"), "\\1", f)
    tbl <- read.delim(file.path(results_dir, f), stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      analysis = label,
      exclusion_type = "within_celltype_skipped",
      n_excluded = nrow(tbl),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(
      analysis = character(0),
      exclusion_type = character(0),
      n_excluded = integer(0)
    ))
  }
  do.call(rbind, rows)
}

assert_seurat_metadata_column <- function(seurat_obj, column, context = "metadata") {
  if (!column %in% colnames(seurat_obj@meta.data)) {
    stop(
      "Column '", column, "' was not found in Seurat ", context, ". ",
      "Check cell_type_label_go_module in project_parameters.Config.yaml. ",
      "Available columns include: ",
      paste(utils::head(colnames(seurat_obj@meta.data), 40L), collapse = ", "),
      if (ncol(seurat_obj@meta.data) > 40L) ", ..." else ""
    )
  }
  invisible(column)
}

############################################################################################################
