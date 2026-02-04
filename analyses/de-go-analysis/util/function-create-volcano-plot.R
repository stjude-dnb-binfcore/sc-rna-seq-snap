############################################################################################################
#' Run Volcano plots function
#'
#' @param all.markers Dataframe of all.markers
#' @param n_value Annotation database (e.g., org.Hs.eg.db)
#' @param plots_dir Directory to save plots
#'
#' @return Volcano plots
#' @export
#'
#' @examples
#' 
crete_volcano_plot <- function(all.markers,
                               n_value,
                               plots_dir) {

  # ---- Optional: per-cluster volcano plots ----
  group_col <- "cluster"        # change to your grouping column if different
  stopifnot(group_col %in% names(all.markers))
  
  by_grp <- split(all.markers, all.markers[[group_col]])
  by_grp <- by_grp[order(names(by_grp))] #Sort alphabetically by group name
  
  for (cl in names(by_grp)) {
    df <- by_grp[[cl]]
    
    # ---- Make sure we have gene symbols as labels (not numeric rownames) ----
    # Prefer a 'gene' column. If it exists, use it; otherwise try to derive from rownames
    if ("gene" %in% names(df)) {
      df$label <- as.character(df$gene)
    } else if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
      # If rownames exist but might be numbers, still coerce to character
      df$label <- as.character(rownames(df))
    } else {
      warning("Skipping ", cl, ": no labels available.")
      next
    }
    
    # ---- Clamp p-values to avoid Inf ----
    eps <- 1e-300
    df$p_val_adj <- pmax(df$p_val_adj, eps)
    
    # ---- Pick top up/down by adj p and effect size ----
    up <- df %>%
      dplyr::filter(avg_log2FC >= 0.25) %>% #dplyr::filter(avg_log2FC > 0.25) %>% dplyr::arrange(p_val_adj)
      dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
      slice_head(n = n_value)
    
    dn <- df %>%
      dplyr::filter(avg_log2FC <= -0.25) %>% # dplyr::filter(avg_log2FC < 0.25) 
      dplyr::arrange(p_val_adj, avg_log2FC) %>%
      slice_head(n = n_value)
    
    # Use gene labels for selection as well
    select_labels_grp <- unique(c(up$label, dn$label))
    
    # Optionally drop NAs
    select_labels_grp <- select_labels_grp[!is.na(select_labels_grp) & nzchar(select_labels_grp)]
    
    cat("Volcano Plot:", cl, "\n")
    p_grp <- print(EnhancedVolcano::EnhancedVolcano(
      df,
      lab            = df$label,            # <- gene names, not numbers
      x              = "avg_log2FC",
      y              = "p_val_adj",
      selectLab      = select_labels_grp,   # <- same label space as 'lab'
      drawConnectors = TRUE,
      parseLabels    = FALSE,     # If labels still parse as math, disable parsing:
      col            = c("black","black","black","red3"),     # Make significant points red only; others black
      pCutoff        = 5e-2,
      FCcutoff       = 0.25,
      title          = paste("Volcano Plot:", cl),
      subtitle       = NULL,
      ylab           = bquote(~-Log[10]~ 'P adj'),
      # Optional label/point tweaks
      labSize        = 3.3,
      pointSize      = 1.2,
      max.overlaps   = 25
    ))
    
    out_png_cl <- file.path(plots_dir, glue::glue("Volcano-plot-{cl}.png"))
    ggplot2::ggsave(filename = out_png_cl, plot = p_grp, width = 8, height = 6, dpi = 300)
  }

}
############################################################################################################
