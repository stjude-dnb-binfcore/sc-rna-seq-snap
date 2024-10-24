###############################################################
#' Function to shuffle data points for plotting
#'
#' @param df 
#' 
#' @return
#' @export
#'
#' @examples
#' 
shuf <- function(df) {
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}


###############################################################
#' Create plot for UMAP for color gradient
#'
#' @param df 
#' @param umap_val 
#' @param color_value
#' @param palette
#' @param title_name
#' 
#' @return
#' @export
#'
#' @examples
create_UMAP_color_gradient <- function(df, umap_val, color_value, palette, title_name) {
  
  # Define color
  low_color_df <- palette %>%
    filter(color_names == "gradient_3")
  high_color_df <- palette %>%
    filter(color_names == "gradient_8")
  na_color_df <- palette %>%
    filter(color_names == "na_color")

  p <- print(ggplot(shuf(df), aes_string(paste0(umap_val, "_1"), paste0(umap_val, "_2"), color = color_value)) + 
               geom_point(size = 0.1) +
               ggtitle(paste(title_name)) +
               xlab("UMAP_1") +
               ylab("UMAP_2") +
               theme_Publication(base_size = 11) + 
               scale_colour_gradient(low = low_color_df$hex_codes, 
                                     high = high_color_df$hex_codes, 
                                     na.value = na_color_df$hex_codes) +
               theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
  return(p)
}


###############################################################
#' Create plot for UMAP for orig.ident
#'
#' @param df 
#' @param umap_val 
#' @param color_value 
#' @param title_name
#'
#' @return
#' @export
#'
#' @examples
create_UMAP_orig_ident <- function(df, umap_val, color_value, title_name) {
  
  # Define color
  colourCount = length(unique(df$orig.ident))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))

  p <- print(ggplot(shuf(df), aes_string(paste0(umap_val, "_1"), paste0(umap_val, "_2"), color = color_value)) + 
               geom_point(size = 0.1) +
               ggtitle(paste(title_name)) +
               xlab("UMAP_1") +
               ylab("UMAP_2") +
               theme_Publication(base_size = 11) + 
               scale_colour_manual(values=getPalette(colourCount)) +
               guides(color=guide_legend(override.aes=list(size=5))) +
               theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  return(p)
}

###############################################################
#' Create plot for UMAP for ID
#'
#' @param df 
#' @param umap_val 
#' @param color_value 
#' @param title_name
#'
#' @return
#' @export
#'
#' @examples
create_UMAP_ID <- function(df, umap_val, color_value, title_name) {
  
  # Define color
  colourCount = length(unique(df$ID))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))
    
  p <- print(ggplot(shuf(df), aes_string(paste0(umap_val, "_1"), paste0(umap_val, "_2"), color = color_value)) + 
               geom_point(size = 0.1) +
               ggtitle(paste(title_name)) +
               xlab("UMAP_1") +
               ylab("UMAP_2") +
               theme_Publication(base_size = 11) + 
               scale_colour_manual(values=getPalette(colourCount)) +
               guides(color=guide_legend(override.aes=list(size=5))) +
               theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  return(p)
}

###############################################################
#' Create plot for UMAP for condition
#'
#' @param df 
#' @param umap_val 
#' @param color_value
#'
#' @return
#' @export
#'
#' @examples
create_UMAP_condition <- function(df, umap_val, color_value) {
  
  # Define color
  # If we are passing unquoted argument for column names, then convert to string with deparse/substitute and use [[ instead of $.
  colourCount = length(unique(df[[color_value]]))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))
  
  p <- print(ggplot(shuf(df), aes_string(paste0(umap_val, "_1"), paste0(umap_val, "_2"), color = color_value)) + 
               geom_point(size = 0.1) +
               ggtitle(paste(color_value)) +
               xlab("UMAP_1") +
               ylab("UMAP_2") +
               theme_Publication(base_size = 11) + 
               scale_colour_manual(values=getPalette(colourCount)) +
               guides(color=guide_legend(override.aes=list(size=5))) +
               theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
  return(p)
}


###############################################################
#' Create plot for UMAP for condition_split
#'
#' @param df 
#' @param umap_val 
#' @param color_value
#'
#' @return
#' @export
#'
#' @examples
create_UMAP_condition_split <- function(df, umap_val, color_value) {
  
  # Define color
  # If we are passing unquoted argument for column names, then convert to string with deparse/substitute and use [[ instead of $.
  colourCount = length(unique(df[[color_value]]))
  getPalette = colorRampPalette(brewer.pal(colourCount, "Dark2"))

  p <- print(ggplot(shuf(df), aes_string(paste0(umap_val, "_1"), paste0(umap_val, "_2"), color = color_value)) + 
               geom_point(size = 0.1) +
               ggtitle(paste(color_value)) +
               xlab("UMAP_1") +
               ylab("UMAP_2") +
               theme_Publication(base_size = 11) + 
               scale_colour_manual(values=getPalette(colourCount)) +
               guides(color=guide_legend(override.aes=list(size=5))) +
               facet_wrap(as.formula(paste("~", color_value))) +
               theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
   return(p)
}


