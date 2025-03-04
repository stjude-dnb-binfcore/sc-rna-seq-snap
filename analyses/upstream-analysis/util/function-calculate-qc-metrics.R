####################################################################################################
#' Function to calculate the lower and upper thresholds for filtering outliers based on median absolute deviation.
#' Takes input of a numeric column from Seurat metadata
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
Find_Outlier_Thershold <- function(metadata) {
  
  # Low threshold 3 mads below the median of the numeric metadata column
  low_outlier_thershold <- median(metadata) - 3 * mad(metadata)
  
  # High threshold 3 mads above the median of the numeric metadata column
  high_outlier_thershold <- median(metadata) + 3 * mad(metadata)
  
  # Return list of low and high thresholds
  return(list(low_outlier_thershold, high_outlier_thershold))
}


####################################################################################################
#' Function to create violin plot
#'
#' @param seurat_obj 
#' @param grouping 
#' @param min_genes
#' @param min_count
#' @param title_name
#' @param palette 
#'
#' @return
#' @export
#'
#' @examples
create_vln_plot <- function(seurat_obj, grouping, min_genes, min_count, title_name, palette) {
  
  # Define color
  vln_plot_color_df <- palette %>%
    filter(color_names == "vln_plot_color")
  min_color_df <- palette %>%
    filter(color_names == "min_color")

  # Plot number of genes 
  genes_num <- plot_distribution(seurat_obj, 
                                 features = "nFeature_RNA",
                                 grouping = grouping) +
    theme_Publication(base_size = 11) +
    xlab("") + 
    ylab("") + 
    geom_hline(yintercept = min_genes, color = min_color_df$hex_codes, size = 1) +
    scale_fill_manual(values = alpha(c(vln_plot_color_df$hex_codes), .4)) +
    annotate(geom = "text", size = 2.5, x = 1, y = min_genes - 250, label = glue::glue("Min genes: {min_genes}"), color = min_color_df$hex_codes)

  # plot number of UMIs
  umi_num <- plot_distribution(seurat_obj, 
                               features = "nCount_RNA",
                               grouping = grouping) +
    theme_Publication(base_size = 11) +
    xlab("") + 
    ylab("") + 
    geom_hline(yintercept = min_count, color = min_color_df$hex_codes, size = 1) +
    scale_fill_manual(values = alpha(c(vln_plot_color_df$hex_codes), .4)) +
    annotate(geom = "text", size = 2.5, x = 1, y = c(-(min_count+100)), label = glue::glue("Min count: {min_count}"), color = min_color_df$hex_codes)

  # plot percent mitochondrial reads
  mito_num <- plot_distribution(seurat_obj, 
                                features = "percent.mito",
                                grouping = grouping) +
    theme_Publication(base_size = 11) + 
    xlab("") + 
    ylab("") + 
    scale_fill_manual(values = alpha(c(vln_plot_color_df$hex_codes), .4)) +
    annotate(geom = "text", size = 2.5, x = 1, y = c(-(mtDNA_pct_default+100)), label = glue::glue("Min mtDNA: {mtDNA_pct_default}"), color = min_color_df$hex_codes)
  
  # get the legend for one of the plots to use as legend for the combined plot
  legend_grid <- get_legend(mito_num)
  
  # Now add the title
  title_name_grid <- ggdraw() +
    draw_label(title_name, fontface = "bold", x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  
  # Combine plots 
  plot_row <- plot_grid(genes_num + theme(legend.position = "none"),
                        umi_num + theme(legend.position = "none"),
                        mito_num + theme(legend.position = "none"),
                        legend_grid, 
                        ncol = 4)
  
  plot <- plot_grid(title_name_grid, plot_row, ncol = 1, rel_heights = c(0.1, 1))
  
  return(plot)
}

###############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#' 
#'
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_1 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open path/name to write figures to 
  # pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname)

  # Make density plots of UMI count, gene count split and colored by sample ID
  plot1 <- ggplot() + geom_density(data = Seurat_obj@meta.data, mapping = aes(x = nCount_RNA, fill = ID), alpha = 0.2) + 
    scale_x_log10() + 
    theme_Publication() + 
    #scale_fill_manual(values = getPalette(colourCount)) +
    ggtitle("nCount_RNA vs Cell Density") + 
    ylab("Cell Density") + 
    theme(aspect.ratio = 1)

  plot2 <- ggplot() + geom_density(data = Seurat_obj@meta.data, mapping = aes(x = nFeature_RNA, fill = ID), alpha = 0.2) +
    scale_x_log10() +
    theme_Publication() + 
    #scale_fill_manual(values = getPalette(colourCount)) +
    ggtitle("nFeature_RNA vs Cell Density") +
    ylab("Cell Density") +
    theme(aspect.ratio = 1)
  
  # get the legend for one of the plots to use as legend for the combined plot
  legend_grid <- get_legend(plot2)
  
  # Combine plots 
  plot_row <- plot_grid(plot1 + theme(legend.position = "none"),
                        plot2 + theme(legend.position = "none"),
                        legend_grid, 
                        ncol = 3)
  plot_grid(plot_row, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(file = fname, width = 14, height = 8, device = "png")
  #dev.off()
  
  return(plot_row)
}
  

############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#' 
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_2 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open pdf path/name to write figures to 
  #pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname)
  
  # Make scatter plot of UMI count vs gene count colored by mitochondrial UMI expression percent
  
  # Define color
  low_color_df <- palette %>%
    filter(color_names == "gradient_3")
  high_color_df <- palette %>%
    filter(color_names == "gradient_8")
  na_color_df <- palette %>%
    filter(color_names == "na_color")
  
  plot_row <- ggplot() + geom_point(data = Seurat_obj@meta.data, mapping = aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mito)) + 
    #scale_color_gradient(low = "grey90", high = "black") + 
    scale_colour_gradient(low = low_color_df$hex_codes, 
                          high = high_color_df$hex_codes, 
                          na.value = na_color_df$hex_codes) +
    theme_Publication() + 
    ggtitle("nFeature_RNA vs nCount_RNA") + 
    theme(aspect.ratio = 1)
  ggsave(file = fname, width = 6, height = 5, device = "png")
  #dev.off()
  
  return(plot_row)
}


############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#' 
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_3 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open pdf path/name to write figures to 
  #pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname)
  
  # Make density plot of mitochondial percent split and colored by sample ID
  plot1 <- ggplot() + geom_density(data = Seurat_obj@meta.data, mapping = aes(x = percent.mito, fill = ID), alpha = 0.2) +
    scale_x_log10() +
    theme_Publication() + 
    #scale_fill_manual(values = getPalette(colourCount)) +
    ggtitle("percent.mito vs Cell Density") + 
    ylab("Cell Density") + 
    theme(aspect.ratio = 1)
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  plot2 <- ggplot() + geom_density(data = Seurat_obj@meta.data, mapping = aes(x = log10GenesPerUMI, fill = ID), alpha = 0.2) +
    theme_Publication() +
    #scale_fill_manual(values = getPalette(colourCount)) +
    ggtitle("log10GenesPerUMI vs Cell Density") + 
    ylab("Cell Density") + 
    theme(aspect.ratio = 1)
  
  # get the legend for one of the plots to use as legend for the combined plot
  legend_grid <- get_legend(plot2)
  
  # Combine plots 
  plot_row <- plot_grid(plot1 + theme(legend.position = "none"),
                        plot2 + theme(legend.position = "none"),
                        legend_grid, 
                        ncol = 3)
  plot_grid(plot_row, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(file = fname, width = 14, height = 8, device = "png")
  #dev.off()
  
  return(plot_row)
}
  

############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#' 
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_4 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open pdf path/name to write figures to 
  #pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname) 
  
  # Make scatter plots of UMI count vs mitochondrial UMI expression percent colored by sample ID and colored by cell density
  plot1 <- FeatureScatter(object = Seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "ID") + 
    theme_Publication() +
    ggtitle("nCount_RNA vs percent.mito") +
    theme(aspect.ratio = 1)
  
  #set.seed(1234)
  plot2 <- FeatureScatter(object = Seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "ID") + 
    geom_bin_2d(bins = 300) + 
    scale_fill_continuous(type = "viridis") + 
    theme_Publication() + 
    ggtitle("nCount_RNA vs percent.mito") +
    theme(aspect.ratio = 1)
  
  # get the legend for one of the plots to use as legend for the combined plot
  legend_grid <- get_legend(plot2)
  
  # Combine plots 
  plot_row <- plot_grid(plot1 + theme(legend.position = "none"),
                        plot2 + theme(legend.position = "none"),
                        legend_grid, 
                        ncol = 3)
  plot_grid(plot_row, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(file = fname, width = 14, height = 8, device = "png")
  #dev.off()
  
  return(plot_row)
} 
  
  
############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#'
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_5 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open pdf path/name to write figures to 
  # pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname) 
  
  # Make scatter plots of UMI count vs gene count colored by sample ID and colored by cell density
  plot1 <- FeatureScatter(object = Seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "ID") + 
    theme_Publication() +
    ggtitle("nCount_RNA vs nFeature_RNA") +
    theme(aspect.ratio = 1)
  
  #set.seed(1234)
  plot2 <- FeatureScatter(object = Seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "ID") + 
    geom_bin_2d(bins = 300) + 
    scale_fill_continuous(type = "viridis") + 
    theme_Publication() + 
    ggtitle("nCount_RNA vs nFeature_RNA") + 
    theme(aspect.ratio = 1)
  
  # get the legend for one of the plots to use as legend for the combined plot
  legend_grid <- get_legend(plot2)
  
  # Combine plots 
  plot_row <- plot_grid(plot1 + theme(legend.position = "none"),
                        plot2 + theme(legend.position = "none"),
                        legend_grid, 
                        ncol = 3)
  plot_grid(plot_row, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(file = fname, width = 14, height = 8, device = "png")
  #dev.off()
  
  return(plot_row)
} 


############################################################################################################
#' Function to make QC plots a the given Seurat object
#' Project_Path is where sub-directories will be created and figures will be saved
#' @param Seurat_obj
#' @param Project_Path
#' @param Figure_Total
#' @param Analysis_Name
#' @param File_Name
#' 
#' @return
#' @export
#'
#' @examples
#' 
Generate_QC_Plots_6 <- function(Seurat_obj, Project_Path, Figure_Total, Analysis_Name, File_Name, palette) {
  
  # Create directory to save figures to
  dir.create(paste(Project_Path, "/", Analysis_Name, sep = ""))
  
  # Open pdf path/name to write figures to 
  # pdf(file = paste(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".pdf", sep = ""), width = length(unique(Seurat_obj$ID))*12, height = 12)
  fname <- paste0(Project_Path, "/", Analysis_Name, "/", Analysis_Name, "_", Figure_Total, "_", File_Name, ".png", sep = "")
  print(fname) 
  
  # Make bar plot of cell count after filtering, split and colored by sample ID, relative to total cell count after filtering
  plot_row <- ggplot() + geom_bar(data = Seurat_obj@meta.data, mapping = aes(x = ID, fill = ID), width = 0.5, alpha = 0.4) + 
    theme_Publication() + 
    ggtitle("Number of Cells") + 
    geom_label(aes(x = Seurat_obj$ID, y = length(Seurat_obj$ID), label = length(Seurat_obj$ID)))
  ggsave(file = fname, width = 6, height = 5, device = "png")
  #dev.off()
  
  return(plot_row)
}

