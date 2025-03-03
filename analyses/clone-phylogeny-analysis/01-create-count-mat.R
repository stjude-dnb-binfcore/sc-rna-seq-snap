#######################################################################################################
# https://broadinstitute.github.io/KrumlovSingleCellWorkshop2020/data-wrangling-scrnaseq-1.html
# https://www.rdocumentation.org/packages/Seurat/versions/3.1.1/topics/GetAssayData
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(yaml)
  library(tidyverse)
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
# Set up directories and paths to file Inputs/Outputs
root_dir <- yaml$root_dir
metadata_dir <- yaml$metadata_dir
assay <- yaml$assay_annotation_module
analysis_dir <- file.path(root_dir, "analyses", "clone-phylogeny-analysis") 
annotation_results_dir <- file.path(root_dir, "analyses", "cell-types-annotation", "results") 
annotations_dir <- yaml$annotations_dir_rshiny_app
annotations_all_results_dir <- file.path(annotation_results_dir, annotations_dir) 

# Input files
seurat_obj_file <- dir(path = annotations_all_results_dir, pattern =  "seurat_obj.*\\.rds", full.names = TRUE, recursive = TRUE)
print(seurat_obj_file)


# Create module_results_dir
module_results_dir <- 
  file.path(analysis_dir, paste0("results"))
if (!dir.exists(module_results_dir)) {
  dir.create(module_results_dir)
}

# Create step_results_dir
step_results_dir <- 
  file.path(module_results_dir, paste0("01-create-count-mat"))
if (!dir.exists(step_results_dir)) {
  dir.create(step_results_dir)
}


#######################################################
# Read metadata file and define `sample_name`
metadata_file <- file.path(metadata_dir, "project_metadata.tsv") # metadata input file

# Read metadata file and define `sample_name`
project_metadata <- read.csv(metadata_file, sep = "\t", header = TRUE)
sample_name <- unique(as.character(project_metadata$ID))
sample_name <- sort(sample_name, decreasing = FALSE)
print(sample_name)

##########################################################################################################
# Read the Seurat object
seurat_obj <- readRDS(seurat_obj_file)

# Split the Seurat object by the 'orig.ident' in the metadata
seurat_obj_list <- SplitObject(seurat_obj, split.by = "orig.ident")
#seurat_obj_list
##########################################################################################################
# Create list 
#seurat_obj_list <- list()
count_mat_list <- list()


for (i in seq_along(sample_name)){

  # Create results_dir per sample
  results_dir <- file.path(step_results_dir, sample_name[i])
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)}
  
  # Read data per sample
  cat("Processing for sample:", sample_name[i], "\n")
  
  # Read the Seurat object
  #seurat_obj_list[[i]] <- readRDS(seurat_obj_file[i]) # Use double brackets to assign to the list
  # Save the Seurat object to the results directory
  #saveRDS(seurat_obj_list[[i]], file = file.path(results_dir, paste0(sample_name[i], "_seurat_obj.rds")))
  #cat("Saved Seurat object for sample:", sample_name[i], "\n")
  
  # Extract count matrix
  count_mat_list[[i]] <- as.matrix(seurat_obj_list[[i]]@assays[[assay]]@counts)

  # Save the count matrix
  saveRDS(count_mat_list[[i]], file = paste0(results_dir, "/", "count_mat.rds"))
  
  cat("Number of cells in count_mat_list[[i]]:", ncol(count_mat_list[[i]]), "\n")
  cat("Number of genes in count_mat_list[[i]]:", nrow(count_mat_list[[i]]), "\n")
  
  cat("Processing complete for sample:", sample_name[i], "\n")
  
}
##########################################################################################################

