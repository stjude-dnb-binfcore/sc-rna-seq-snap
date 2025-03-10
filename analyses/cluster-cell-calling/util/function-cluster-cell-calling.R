############################################################################################################
#' Function to calculate cluster averages
#' @param metadata
#' @param data
#' @param resolution_list
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
#' 
calculate_cluster_average <- function(metadata, data, resolution_list, assay) {
  
  use_resolution = glue::glue("{assay}_snn_res.{resolution_list}")
  
  # get relevant metadata
  metadata <- metadata %>%
    select("cell", use_resolution)
  
  # manipulate data to merge with metadata
  # genes are already columns, transpose, so genes are rows
  # remove all genes with less than 50 cells
  xx <- tabulate(data@i + 1)
  rows_to_keep <- rownames(data)[xx > 49]
  data <- data[rows_to_keep,]
  
  # split up sparse matrix by clusters and then calculate average of each df individually
  current_mean <- data.frame(rownames(data))
  colnames(current_mean) <- "gene"
  
  for (l in 1:length(unique(metadata[,eval(use_resolution)]))) {
    l <- l-1
    #determine cells belonging to each cluster
    cells <- subset(metadata, metadata[,eval(use_resolution)] == l)[,1]
    # subset for each cluster
    data_pre <- as.data.frame(data[,cells])
    current_mean$cluster <- Matrix::rowMeans(data_pre)
    colnames(current_mean)[2+l] <- paste0("cluster_", l)
  }
  
  return(current_mean)
}


############################################################################################################
#' Function to cluster the cells based on presence/absence of a `resolution_list`
#' @param seurat_obj
#' @param reduction_value
#' @param num_dim
#' @param plots_dir
#' @param resolution_list
#' @param resolution_list_default
#' @param algorithm_value
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
#' 
cluster_cell_calling <- function(seurat_obj, reduction_value, num_dim, plots_dir, resolution_list, resolution_list_default, algorithm_value, assay) {
  
  set.seed(1234) # Make code reproducible

  # If no resolution_list is provided in function call, then use list below (0...10)
  # Generate clustering for each resolution amd generate a figure
  if (is.null(resolution_list)) {
    for (res in resolution_list_default) {
      message("FindClusters is being calculated for resolution: ", res, "\n")
      seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_value, dims = 1:num_dim, verbose = TRUE)
      seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = algorithm_value, method = "igraph", random.seed = 0)
      message("DimPlot is being created for resolution: ", res, "\n")
      
      name <- paste0(plots_dir, "/", "Cluster_resolutions_", res, ".png", sep = "")
      print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution: ", res)))
      ggsave(file = name, width = 12, height = 6, device = "png")
      }

    # Generate cluster tree for all resolutions (to aid in finding "stable" resolution)
    cat("Clustree is being created\n")
    name <- paste0(plots_dir, "/", "Cluster_tree.png", sep = "")
    print(clustree(seurat_obj, prefix = glue::glue("{assay}_snn_res.")))
    ggsave(file = name, width = 25, height = 20, device = "png")
    
    name <- paste0(plots_dir, "/", "Cluster_tree_sc3_stability.png", sep = "")
    print(clustree(seurat_obj, prefix = glue::glue("{assay}_snn_res."), node_colour = "sc3_stability"))
    ggsave(file = name, width = 25, height = 20, device = "png")

    print_message <- print("we use multiple resolutions by default for now")
    
    #### #### ####
    #metadata <- as_data_frame_seurat(seurat_obj, metadata = TRUE)
  
    # Identify columns with a '.1' suffix
    cols_to_remove <- grep("\\.1$", colnames(seurat_obj@meta.data), value = TRUE)
    
    # Exclude columns that match the specific patterns (e.g., {assay}_snn_res.0.1, {assay}_snn_res.1, {assay}_snn_res.10)
    cols_to_remove <- cols_to_remove[!grepl(glue::glue("^{assay}_snn_res\\.0\\.1$"), cols_to_remove) & 
                                       !grepl(glue::glue("^{assay}_snn_res\\.1$"), cols_to_remove) &
                                       !grepl(glue::glue("^{assay}_snn_res\\.10$"), cols_to_remove)]
    
    # Remove the columns
    seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% cols_to_remove]
    
    metadata <- seurat_obj@meta.data
    write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_integrated_{integration_method}_clusters_all.tsv")))
    #### #### ####
    
    saveRDS(seurat_obj, file = paste0(results_dir, "/", glue::glue("seurat_obj_integrated_{integration_method}_clusters_all.rds"))) # save object
    
    # If resolution_list is provided, only generate clustering and figures for those resolutions
    } else if (length(resolution_list) > 1) {
      for (res in resolution_list) {
        message("FindClusters is being calculated for resolution: ", res, "\n")
        seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_value, dims = 1:num_dim, verbose = TRUE)
        seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = algorithm_value, method = "igraph", random.seed = 0)
        message("DimPlot is being created for resolution: ", res, "\n")
        name <- paste0(plots_dir, "/", "Cluster_resolutions_", res, ".png", sep = "")
        print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution: ", res)))
        ggsave(file = name, width = 12, height = 6, device = "png")
      }
      #dev.off()

      
      # Generate cluster tree for all resolutions (to aid in finding "stable" resolution)
      cat("Clustree is being created\n")
      #pdf(file = paste0(plots_dir, "/", "Cluster_tree.pdf", sep = ""), width = 25, height = 20) #60
      name <- paste0(plots_dir, "/", "Cluster_tree.png", sep = "")
      print(clustree(seurat_obj, prefix = glue::glue("{assay}_snn_res.")))
      ggsave(file = name, width = 12, height = 8, device = "png")
      
      name <- paste0(plots_dir, "/", "Cluster_tree_sc3_stability.png", sep = "")
      print(clustree(seurat_obj, prefix = glue::glue("{assay}_snn_res."), node_colour = "sc3_stability"))
      ggsave(file = name, width = 12, height = 8, device = "png")
      #dev.off()
      
      print_message <- print("we use multiple resolutions for now")
        
      #### #### ####
      #metadata <- as_data_frame_seurat(seurat_obj, metadata = TRUE)
      
      # Identify columns with a '.1' suffix
      cols_to_remove <- grep("\\.1$", colnames(seurat_obj@meta.data), value = TRUE)
      
      ## Exclude columns that match the specific patterns (e.g., {assay}_snn_res.0.1, {assay}_snn_res.1, {assay}_snn_res.10)
      cols_to_remove <- cols_to_remove[!grepl(glue::glue("^{assay}_snn_res\\.0\\.1$"), cols_to_remove) & 
                                         !grepl(glue::glue("^{assay}_snn_res\\.1$"), cols_to_remove) &
                                         !grepl(glue::glue("^{assay}_snn_res\\.10$"), cols_to_remove)]
      
      # Remove the columns
      seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% cols_to_remove]
      
      metadata <- seurat_obj@meta.data
      write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_integrated_{integration_method}_clusters_all.tsv")))
      #### #### ####
        
        saveRDS(seurat_obj, file = paste0(results_dir, "/", glue::glue("seurat_obj_integrated_{integration_method}_clusters_all.rds"))) # save object
        
      # If a single resolution is provided, only use that resolution to generate clustering and figure
      # Don't generate cluster tree
      } else if (length(resolution_list) == 1) {
        for (res in resolution_list) {
          message("FindClusters is being calculated for resolution: ", res, "\n")
          seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_value, dims = 1:num_dim, verbose = TRUE)
          seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = algorithm_value, method = "igraph", random.seed = 0)
          message("DimPlot is being created for resolution: ", res, "\n")
          name <- paste0(plots_dir, "/", "Cluster_resolutions_", res, ".png", sep = "")
          print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution: ", res)))
          ggsave(file = name, width = 12, height = 6, device = "png")
        }
        #dev.off()
        
        print_message <- print("we use the sigle resolution that fits the data best")
        # Calculate average
        counts_matrix <- seurat_obj@assays$assay@counts 
        #metadata <- as_data_frame_seurat(seurat_obj, metadata = TRUE)
        # Identify columns with a '.1' suffix
        cols_to_remove <- grep("\\.1$", colnames(seurat_obj@meta.data), value = TRUE)
        
        # Exclude columns that match the specific patterns (e.g., {assay}_snn_res.0.1, {assay}_snn_res.1, {assay}_snn_res.10)
        cols_to_remove <- cols_to_remove[!grepl(glue::glue("^{assay}_snn_res\\.0\\.1$"), cols_to_remove) & 
                                           !grepl(glue::glue("^{assay}_snn_res\\.1$"), cols_to_remove) &
                                           !grepl(glue::glue("^{assay}_snn_res\\.10$"), cols_to_remove)]
        
        # Remove the columns
        seurat_obj@meta.data <- seurat_obj@meta.data[, !colnames(seurat_obj@meta.data) %in% cols_to_remove]
        
        metadata <- seurat_obj@meta.data
        avg_res <- calculate_cluster_average(metadata = metadata, 
                                             data = counts_matrix,
                                             resolution_list = resolution_list)
        write_tsv(avg_res, file = paste0(results_dir, "/", glue::glue("{assay}_snn_res.{resolution_list}-avg.tsv"))) 
        #### #### ####
        write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_integrated_{integration_method}_clusters_{resolution_list}.tsv")))
        #### #### ####
        
        saveRDS(seurat_obj, file = paste0(results_dir, "/", glue::glue("seurat_obj_integrated_{integration_method}_clusters_{resolution_list}.rds"))) # save object
        }
  
  return(seurat_obj) }

############################################################################################################
