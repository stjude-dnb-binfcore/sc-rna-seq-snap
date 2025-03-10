#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(Seurat)
  library(scooter)
  })

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)
#################################################################################
# Parameters
root_dir <- yaml$root_dir
PROJECT_NAME <- yaml$PROJECT_NAME
method <- yaml$method

# Set up directories and paths to root_dir and analysis_dir
analysis_dir <- file.path(root_dir, "analyses", "cell-types-annotation") 
module_results_dir <- file.path(analysis_dir, "results")

broad_SingleR_results_dir <- file.path(module_results_dir, "01_cell_types_annotation_SingleR_broad") 
fine_SingleR_results_dir <- file.path(module_results_dir, "02_cell_types_annotation_SingleR_fine") 
gene_markers_results_dir <- file.path(module_results_dir, "03_cell_types_annotation_gene_markers") 

# Create dir
results_dir <- file.path(module_results_dir, "04_cell_types_annotations_all")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}

################################################################################################################
### Generate final object ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################################################################################################

if (method == "all"){
  # Input files
  broad_SingleR_file <- file.path(broad_SingleR_results_dir, "seurat_obj_SingleR_broad.rds")
  fine_SingleR_file <- file.path(fine_SingleR_results_dir, "metadata.tsv")
  gene_markers_file <- file.path(gene_markers_results_dir, "metadata.tsv")

  # Read seurat object #########################
  ##############################################
  seurat_obj <- readRDS(broad_SingleR_file)
  
  # Rename the column in the metadata of the Seurat object
  seurat_obj@meta.data$pruned.labels.broad <- seurat_obj@meta.data$pruned.labels
  
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    rename_with(~ paste("broad_", ., sep = ""), matches("^scores\\."))

  # We want to attach the metadata related to the cell type annotation from other methods to the ones of the object
  # We randomly chose the object from the first method, `SingleR_broad`

  # Read metadata `fine_SingleR`
 fine_SingleR_df <- readr::read_tsv(fine_SingleR_file, guess_max = 100000, show_col_types = FALSE) %>%
   mutate(pruned.labels.fine = pruned.labels) %>%
   
   # Select for columns to use for join to the object
   select(cell, pruned.labels.fine, singler.fine, matches("^scores\\.")) %>% # start with
   rename_with(~ paste("fine_", ., sep = ""), matches("^scores\\."))
 
 
  # Read metadata `gene_markers`
  new_metadata <- readr::read_tsv(gene_markers_file, guess_max = 100000, show_col_types = FALSE) %>%
    
    # Select for columns to use for join to the object
    select(cell, predicted.cell.signature.ident, matches("\\.score1$")) %>% #end with
    rename_with(~ paste("gene_markers_", ., sep = ""), matches("\\.score1$")) %>%
  
    # Join df
    left_join(fine_SingleR_df) %>%
    select(!cell)

    # Add metadata
    seurat_obj <- AddMetaData(seurat_obj, metadata = new_metadata)
    
    } else if (method == "singler"){
      
      # Input files
      broad_SingleR_file <- file.path(broad_SingleR_results_dir, "seurat_obj_SingleR_broad.rds")
      fine_SingleR_file <- file.path(fine_SingleR_results_dir, "metadata.tsv")

      # Read seurat object #########################
      ##############################################
      seurat_obj <- readRDS(broad_SingleR_file)
      
      # Rename the column in the metadata of the Seurat object
      seurat_obj@meta.data$pruned.labels.broad <- seurat_obj@meta.data$pruned.labels
      
      seurat_obj@meta.data <- seurat_obj@meta.data %>%
        rename_with(~ paste("broad_", ., sep = ""), matches("^scores\\."))
      
      # We want to attach the metadata related to the cell type annotation from other methods to the ones of the object
      # We randomly chose the object from the first method, `SingleR_broad`
      
      # Read metadata `fine_SingleR`
      new_metadata <- readr::read_tsv(fine_SingleR_file, guess_max = 100000, show_col_types = FALSE) %>%
        mutate(pruned.labels.fine = pruned.labels) %>%
        
        # Select for columns to use for join to the object
        select(pruned.labels.fine, singler.fine, matches("^scores\\.")) %>% # start with
        rename_with(~ paste("fine_", ., sep = ""), matches("^scores\\."))
      
      # Check if the number of rows in new_metadata matches the number of cells in the Seurat object
      if (nrow(new_metadata) == ncol(seurat_obj)) {
        # Add metadata to Seurat object
        seurat_obj <- AddMetaData(seurat_obj, metadata = new_metadata)

      } else {
        # If the number of rows doesn't match, print a warning
        warning("Number of rows in new_metadata doesn't match the number of cells in seurat_obj. Label metadata will not be added to the Seurat object.")
      }
    }
  
  
#############################################################################
# Save output files #########################
reduction_names <- c(paste0("umap")) # Export the reductions to Seurat
metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)

write_tsv(metadata, file = paste0(results_dir, "/", "metadata", ".tsv")) # Save metadata
saveRDS(seurat_obj, file = paste0(results_dir, "/", "seurat_obj_cell_types_annotations_all.rds"))
################################################################################################################   
