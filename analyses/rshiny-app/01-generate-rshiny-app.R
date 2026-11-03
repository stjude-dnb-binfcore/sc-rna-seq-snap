#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(reticulate)
  library(gridExtra)
  library(glue)
  library(R.utils)
  library(shiny)
  library(shinyhelper)
  library(DT)
  library(magrittr)
  library(ggdendro)
  library(ShinyCell)
  library(RColorBrewer)
})

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file (UTF-8; same approach as de-go-analysis / project-updates)
config_lines <- readLines(configFile, encoding = "UTF-8", warn = FALSE)
yaml <- read_yaml(text = paste(config_lines, collapse = "\n"))
#################################################################################
# Parameters
root_dir <- yaml$root_dir
PROJECT_NAME <- yaml$PROJECT_NAME
PI_NAME <- yaml$PI_NAME
condition_value <- yaml$condition_value
assay <- yaml$assay_filter_object
annotations_dir <- yaml$annotations_dir_rshiny_app
annotations_filename <- yaml$annotations_filename_rshiny_app
reduction_value <- yaml$reduction_value_annotation_module
resolution_list_dego <- yaml$resolution_list_find_markers
dego_mode_rshiny <- yaml$dego_condition_comparison_mode
if (is.null(dego_mode_rshiny) || !nzchar(as.character(dego_mode_rshiny))) {
  dego_mode_rshiny <- "global_interaction"
}

is_dego_shiny_yes <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) {
    return(FALSE)
  }
  if (isTRUE(x)) {
    return(TRUE)
  }
  identical(toupper(trimws(as.character(x))), "YES")
}
include_de_go_rshiny <- is_dego_shiny_yes(yaml$include_de_go_analysis_rshiny_app_module)


# Set up directories and paths to root_dir and analysis_dir
analysis_dir <- file.path(root_dir, "analyses") 
module_dir <- file.path(analysis_dir, "rshiny-app") 

annotation_results_dir <- file.path(analysis_dir, "cell-types-annotation", "results") 
#broad_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "01_annotations_broad") 
#fine_SingleR_results_dir <- file.path(annotation_results_dir, "01_cell_types_annotation_SingleR", "02_annotations_fine") 
annotations_all_results_dir <- file.path(annotation_results_dir, annotations_dir) 

# Input files
#broad_SingleR_file <- file.path(broad_SingleR_results_dir, "seurat_obj_SingleR_broad.rds")
#fine_SingleR_file <- file.path(fine_SingleR_results_dir, "seurat_obj_SingleR_fine.rds")
annotations_all_file <- file.path(annotations_all_results_dir, annotations_filename)

# Create results_dir
results_dir <- file.path(module_dir, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}


################################################################################
### library(ShinyCell) ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
################################################################################
# It was downloaded from here: https://github.com/SGDDNB/ShinyCell/blob/master/R
# We modified the function to account for any type of assay
source(paste0(module_dir, "/util/makeShinyFiles_assay.R"))
source(paste0(module_dir, "/util/function-attach-dego-shiny-app.R"))

################################################################################################################
### Generate R shiny app ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
################################################################################################################
#cat("Beginning to process results from", "broad_SingleR_file", "\n")
#seu1 <- readRDS(broad_SingleR_file)
#scConf1 <- createConfig(seu1)
#makeShinyFiles_assay(seu1, scConf1, shiny.prefix = "sc1", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))


#cat("Beginning to process results from", "fine_SingleR_file", "\n")
#seu2 <- readRDS(fine_SingleR_file)
#scConf2 <- createConfig(seu2)
#makeShinyFiles_assay(seu2, scConf2, shiny.prefix = "sc2", shiny.dir = paste(results_dir, "shinyApp", sep = "/"))

cat("Beginning to process results from", "annotations_file", "\n")
seu1 <- readRDS(annotations_all_file)


# Check if the specified dimensionality reduction exists in the Seurat object
if (reduction_value %in% names(seu1@reductions)) {
  
  # Retrieve the key prefix for the specified reduction (e.g., "PC_", "UMAP_")
  reduc_key <- seu1[[reduction_value]]@key
  message("Setting default reduction to ", reduction_value)
  
  # Construct names for dimension 1 and 2 (e.g., "PC_1", "PC_2")
  dim_1 <- paste0(reduc_key, "1")
  dim_2 <- paste0(reduc_key, "2")
  
  # Identify all metadata columns associated with the specified reduction
  cols_to_remove <- grep(reduc_key, colnames(seu1@meta.data), value = TRUE)
  
  # Remove the identified columns from the metadata
  seu1@meta.data <- seu1@meta.data[, !colnames(seu1@meta.data) %in% cols_to_remove]
  
} else {
  # If the reduction is not present, notify that the default will be set automatically
  message(reduction_value, " not found in Seurat object, default reduction will be determined by the function")
}

scConf1 <- createConfig(seu1, 
                        meta.to.include = NA, # Include all metadata (or specify if you want a subset)
                        maxLevels = 150)     # Use the number of unique levels

makeShinyFiles_assay(seu1, 
                     scConf1, 
                     shiny.prefix = "sc1", 
                     #default.dimred = c("umap1", "umap2"),
                     default.dimred = c(dim_1, dim_2),
                     shiny.dir = paste(results_dir, "shinyApp", sep = "/"))

cat("Make R shiny app for all files", "\n")
shiny_app_dir <- file.path(results_dir, "shinyApp")
makeShinyCodesMulti(
  shiny.title = PROJECT_NAME,
  shiny.footnotes = PI_NAME,
  shiny.prefix = c("sc1"),
  shiny.headers = c("Annotations"),
  #shiny.prefix = c("sc1", "sc2"),
  #shiny.headers = c("broad_SingleR", "fine_SingleR"),
  shiny.dir = shiny_app_dir)

if (include_de_go_rshiny) {
  attach_de_go_reports_to_shiny_app(
    shiny_dir = shiny_app_dir,
    root_dir = root_dir,
    dego_mode = dego_mode_rshiny,
    resolution_list = resolution_list_dego
  )
} else {
  message(
    "de-go reports not added to Shiny app ",
    "(set include_de_go_analysis_rshiny_app_module: \"YES\" in project_parameters.Config.yaml)"
  )
}
################################################################################################################
