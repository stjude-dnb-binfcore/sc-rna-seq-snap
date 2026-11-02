############################################################################################################
#' Original de-go condition DE: cell_type x condition identity vs all other cells
#'
#' @export
run_global_interaction_de_for_column <- function(seurat_obj,
                                                 cell_type_col,
                                                 condition_col,
                                                 output_suffix,
                                                 results_dir,
                                                 resolution_list,
                                                 assay,
                                                 min_cells_per_group = 25,
                                                 logfc.threshold = 0.25,
                                                 min.pct = 0.1,
                                                 test.use = "wilcox",
                                                 return.thresh = 0.01) {
  Idents(seurat_obj) <- interaction(
    seurat_obj@meta.data[[cell_type_col]],
    seurat_obj@meta.data[[condition_col]],
    drop = TRUE
  )
  message(
    "FindAllMarkers (global interaction): column=", condition_col,
    " (", output_suffix, "); ", length(levels(Idents(seurat_obj))), " identities"
  )
  all.markers <- FindAllMarkers(
    seurat_obj,
    assay = assay,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    min.pct = min.pct,
    return.thresh = return.thresh
  )
  res <- finalize_global_interaction_markers(
    markers = all.markers,
    idents_counts = table(Idents(seurat_obj)),
    min_cells_per_group = min_cells_per_group,
    results_dir = results_dir,
    resolution_list = resolution_list,
    output_suffix = output_suffix
  )
  invisible(res$path)
}

normalize_dego_comparison_mode <- function(mode) {
  if (is.null(mode) || mode %in% c("NULL", ".", "")) {
    return("global_interaction")
  }
  mode <- tolower(as.character(mode))
  mode <- gsub("-", "_", mode)
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

run_dego_mode_needs_precompute <- function(mode) {
  mode %in% c("within_celltype", "both")
}

run_dego_mode_runs_global_in_notebook <- function(mode) {
  identical(mode, "global_interaction")
}

############################################################################################################
