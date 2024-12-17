#' Function to quantify cell type annotations fractions - Broad cell type annotation
#' @param df
#' @param condition
#' @param color_df
#' @param title_value
#' 
#' @return
#' @export
#'
#' @examples
#' 
cell_type_fractions_broad <- function(df, condition, color_df, title_value) {
  
  # Calculate fractions
  count <- df %>% group_by_at(condition) %>%
    mutate(condition=row_number()) %>% 
    dplyr::count(singler.broad)
  
  # Order cell types
  cell_type_order <- unique(as.character(count$singler.broad))
  cell_type_order <- sort(cell_type_order, decreasing = FALSE)
  
  # Define label for plots
  count$singler.broad <- factor(count$singler.broad, levels = cell_type_order)
  
  # Plot
  p <- ggplot(count, aes(x = count[[condition]], y = n, fill = singler.broad)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = color_df) +
    theme_Publication(base_size = 11) + 
    xlab(glue::glue("{condition}")) +
    ylab("Percent Cell Type") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggtitle(glue::glue("{title_value} cell type fractions per {condition}")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(p)
}



############################################################################################################
#' Function to quantify cell type annotations fractions - Fine cell type annotation
#' @param df
#' @param condition
#' @param color_df
#' @param title_value
#' 
#' @return
#' @export
#'
#' @examples
#' 
cell_type_fractions_fine <- function(df, condition, color_df, title_value) {
  
  # Calculate fractions
  count <- df %>% group_by_at(condition) %>%
    mutate(condition=row_number()) %>% 
    dplyr::count(singler.fine)
  
  # Order cell types
  cell_type_order <- unique(as.character(count$singler.fine))
  cell_type_order <- sort(cell_type_order, decreasing = FALSE)
  
  # Define label for plots
  count$singler.fine <- factor(count$singler.fine, levels = cell_type_order)
  
  # Plot
  p <- ggplot(count, aes(x = count[[condition]], y = n, fill = singler.fine)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = color_df) +
    theme_Publication(base_size = 11) + 
    xlab(glue::glue("{condition}")) +
    ylab("Percent Cell Type") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggtitle(glue::glue("{title_value} cell type fractions per {condition}")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(p)
}
############################################################################################################
