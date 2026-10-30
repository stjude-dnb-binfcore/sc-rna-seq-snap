suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

############################################################################################################
#' Differential expression: condition contrast within each cell type
#'
#' Subsets to one cell type at a time and runs FindMarkers between two condition
#' levels (ident.1 vs ident.2). Cell types with fewer than min_cells in either
#' group are skipped.
#'
#' @param seurat_obj Seurat object
#' @param cell_type_col Metadata column for cell type labels
#' @param condition_col Metadata column for experimental condition
#' @param assay Assay name
#' @param ident1 Level in ident.1 (higher expression = up in ident.1)
#' @param ident2 Reference level (ident.2)
#' @param min_cells_per_group Minimum cells required in ident.1 and ident.2 per cell type
#' @param logfc.threshold Passed to FindMarkers
#' @param min.pct Passed to FindMarkers
#' @param test.use Passed to FindMarkers
#' @param return.thresh Passed to FindMarkers
#' @param parallel_workers Number of workers for future_lapply (1 = sequential)
#'
#' @return List with elements markers (data.frame), skipped (data.frame)
#' @export
run_within_celltype_de <- function(seurat_obj,
                                   cell_type_col,
                                   condition_col,
                                   assay,
                                   ident1,
                                   ident2,
                                   min_cells_per_group = 25,
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1,
                                   test.use = "wilcox",
                                   return.thresh = 0.01,
                                   parallel_workers = 1) {
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cell type column not found in metadata: ", cell_type_col)
  }
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("Condition column not found in metadata: ", condition_col)
  }

  md <- seurat_obj@meta.data
  cell_types <- sort(unique(as.character(md[[cell_type_col]])))
  cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]

  run_one_cell_type <- function(ct) {
    cells_use <- rownames(md)[
      md[[cell_type_col]] == ct &
        !is.na(md[[condition_col]]) &
        md[[condition_col]] %in% c(ident1, ident2)
    ]
    n1 <- sum(md[cells_use, condition_col] == ident1, na.rm = TRUE)
    n2 <- sum(md[cells_use, condition_col] == ident2, na.rm = TRUE)

    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      return(list(
        markers = NULL,
        skipped = data.frame(
          cell_type = ct,
          condition_column = condition_col,
          ident1 = ident1,
          ident2 = ident2,
          n_ident1 = n1,
          n_ident2 = n2,
          reason = "below_min_cells_per_group",
          stringsAsFactors = FALSE
        )
      ))
    }

    sub_obj <- subset(seurat_obj, cells = cells_use)
    Idents(sub_obj) <- sub_obj@meta.data[[condition_col]]
    mk <- tryCatch(
      FindMarkers(
        sub_obj,
        assay = assay,
        ident.1 = ident1,
        ident.2 = ident2,
        logfc.threshold = logfc.threshold,
        test.use = test.use,
        min.pct = min.pct,
        return.thresh = return.thresh
      ),
      error = function(e) {
        message("FindMarkers failed for ", ct, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(mk) || nrow(mk) == 0) {
      return(list(
        markers = NULL,
        skipped = data.frame(
          cell_type = ct,
          condition_column = condition_col,
          ident1 = ident1,
          ident2 = ident2,
          n_ident1 = n1,
          n_ident2 = n2,
          reason = "no_markers_or_error",
          stringsAsFactors = FALSE
        )
      ))
    }

    mk$gene <- rownames(mk)
    mk$cluster <- ct
    mk$condition_column <- condition_col
    mk$ident1 <- ident1
    mk$ident2 <- ident2
    mk$n_ident1 <- n1
    mk$n_ident2 <- n2

    list(markers = mk, skipped = NULL)
  }

  if (parallel_workers > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("future.apply not installed; running cell types sequentially")
      parallel_workers <- 1L
    }
  }
  if (parallel_workers > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = parallel_workers)
    res_list <- future.apply::future_lapply(
      cell_types,
      run_one_cell_type,
      future.seed = TRUE
    )
  } else {
    res_list <- lapply(cell_types, run_one_cell_type)
  }

  marker_parts <- lapply(res_list, function(x) x$markers)
  marker_parts <- marker_parts[!vapply(marker_parts, is.null, logical(1L))]
  markers <- if (length(marker_parts) == 0L) {
    tibble::tibble()
  } else {
    bind_rows(marker_parts)
  }

  skipped_parts <- lapply(res_list, function(x) x$skipped)
  skipped_parts <- skipped_parts[!vapply(skipped_parts, is.null, logical(1L))]
  skipped <- if (length(skipped_parts) == 0L) {
    tibble::tibble()
  } else {
    bind_rows(skipped_parts)
  }

  list(markers = markers, skipped = skipped)
}

############################################################################################################
#' Resolve ident.1 / ident.2 for a condition metadata column
#'
#' @param seurat_obj Seurat object
#' @param condition_col Condition metadata column
#' @param ident1 User ident.1 or NULL
#' @param ident2 User ident.2 or NULL
#'
#' @return Named list ident1, ident2, levels
resolve_de_condition_idents <- function(seurat_obj, condition_col, ident1 = NULL, ident2 = NULL) {
  levs <- sort(unique(as.character(seurat_obj@meta.data[[condition_col]])))
  levs <- levs[!is.na(levs) & nzchar(levs)]

  if (length(levs) < 2L) {
    stop(
      "Condition column ", condition_col, " must have at least two levels; found: ",
      paste(levs, collapse = ", ")
    )
  }

  if (!is.null(ident1) && !is.null(ident2)) {
    if (!ident1 %in% levs || !ident2 %in% levs) {
      stop(
        "dego_ident1/dego_ident2 must be levels of ", condition_col, ": ",
        paste(levs, collapse = ", ")
      )
    }
    return(list(ident1 = ident1, ident2 = ident2, levels = levs))
  }

  if (length(levs) == 2L) {
    return(list(ident1 = levs[2L], ident2 = levs[1L], levels = levs))
  }

  stop(
    "Condition column ", condition_col, " has ", length(levs),
    " levels. Set dego_ident1 and dego_ident2 in project_parameters.Config.yaml (",
    paste(levs, collapse = ", "), ")."
  )
}

############################################################################################################
#' Run within-cell-type DE for one condition metadata column and write TSV outputs
#'
#' @param seurat_obj Seurat object
#' @param output_suffix Suffix for output files (e.g. condition_value1)
#' @param ... Passed to run_within_celltype_de and resolve helpers via args list
#' @param results_dir Directory for TSV files
#' @param resolution_list Resolution label used in filenames
run_within_celltype_de_for_column <- function(seurat_obj,
                                              condition_col,
                                              output_suffix,
                                              results_dir,
                                              resolution_list,
                                              cell_type_col,
                                              assay,
                                              ident1 = NULL,
                                              ident2 = NULL,
                                              min_cells_per_group = 25,
                                              parallel_workers = 1,
                                              logfc.threshold = 0.25,
                                              min.pct = 0.1,
                                              test.use = "wilcox",
                                              return.thresh = 0.01) {
  idents <- resolve_de_condition_idents(seurat_obj, condition_col, ident1, ident2)
  message(
    "Within-cell-type DE: column=", condition_col,
    " (", output_suffix, ") ident.1=", idents$ident1,
    " vs ident.2=", idents$ident2
  )

  res <- run_within_celltype_de(
    seurat_obj = seurat_obj,
    cell_type_col = cell_type_col,
    condition_col = condition_col,
    assay = assay,
    ident1 = idents$ident1,
    ident2 = idents$ident2,
    min_cells_per_group = min_cells_per_group,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    test.use = test.use,
    return.thresh = return.thresh,
    parallel_workers = parallel_workers
  )

  out_markers <- file.path(
    results_dir,
    paste0("Res_", resolution_list, "_Markers_within_celltype_", output_suffix, ".tsv")
  )
  out_skipped <- file.path(
    results_dir,
    paste0("Res_", resolution_list, "_Skipped_celltypes_", output_suffix, ".tsv")
  )

  if (nrow(res$markers) > 0) {
    write.table(res$markers, file = out_markers, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    message("No within-cell-type markers written for ", output_suffix)
    write.table(
      data.frame(note = "No cell types passed filters or returned markers"),
      file = out_markers,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (nrow(res$skipped) > 0) {
    write.table(res$skipped, file = out_skipped, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  invisible(list(markers_path = out_markers, skipped_path = out_skipped, result = res))
}

############################################################################################################
