#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
#.libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(rmarkdown)
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

# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "project-updates") 
metadata_dir <- yaml$metadata_dir
metadata_file <- yaml$metadata_file

PROJECT_NAME <- yaml$PROJECT_NAME
genome_reference <- yaml$genome_name
min_genes <- yaml$min_genes
min_count <- yaml$min_count
mtDNA_pct_default <- yaml$mtDNA_pct_default
resolution_clustering_module <- yaml$resolution_clustering_module
resolution_values <- yaml$resolution_list_default_clustering_module
resolution_list_find_markers <- yaml$resolution_list_find_markers
n_value_find_markers <- yaml$n_value_find_markers
method_annotation_module <- yaml$method_annotation_module
annotation_file1 <- yaml$gene_markers_file_name_annotation_module
annotation_file2 <- yaml$reference_file_name_annotation_module
integration_method <-  yaml$integration_method_clustering_module

# Set variables for `project-updates` module
results_filepath <- yaml$results_filepath
mac_results_filepath <- yaml$mac_results_filepath
cohort_value <- yaml$cohort_value

is_yes_param <- function(x) {
  identical(toupper(as.character(x)), "YES")
}

include_de_go_project_updates <- is_yes_param(yaml$include_de_go_analysis_project_updates_module)
dego_config_util <- file.path(root_dir, "analyses", "de-go-analysis", "util", "function-de-go-config.R")
if (include_de_go_project_updates && file.exists(dego_config_util)) {
  source(dego_config_util)
  dego_mode_project_updates <- normalize_dego_mode_label(yaml$dego_condition_comparison_mode)
} else {
  dego_mode_project_updates <- "global_interaction"
}

report_params <- list(
  cellranger_parameters = yaml$cellranger_parameters,
  root_dir = yaml$root_dir,
  metadata_dir = yaml$metadata_dir,
  metadata_file = yaml$metadata_file,
  results_filepath = yaml$results_filepath,
  mac_results_filepath = yaml$mac_results_filepath,
  cohort_value = yaml$cohort_value,
  genome_reference = yaml$genome_name,
  min_genes = yaml$min_genes,
  min_count = yaml$min_count,
  mtDNA_pct_default = yaml$mtDNA_pct_default,
  resolution_clustering_module = yaml$resolution_clustering_module,
  resolution_values = yaml$resolution_list_default_clustering_module,
  resolution_list_find_markers = yaml$resolution_list_find_markers,
  n_value_find_markers = yaml$n_value_find_markers,
  method_annotation_module = yaml$method_annotation_module,
  annotation_file1 = yaml$gene_markers_file_name_annotation_module,
  annotation_file2 = yaml$reference_file_name_annotation_module,
  integration_method = yaml$integration_method_clustering_module,
  include_de_go_analysis_project_updates_module = yaml$include_de_go_analysis_project_updates_module,
  dego_condition_comparison_mode = dego_mode_project_updates,
  cell_type_label_go_module = yaml$cell_type_label_go_module,
  min_cells_per_group_dego_module = yaml$min_cells_per_group_dego_module,
  condition_value1 = yaml$condition_value1,
  condition_value2 = yaml$condition_value2,
  condition_value3 = yaml$condition_value3,
  PROJECT_NAME = yaml$PROJECT_NAME,
  PI_NAME = yaml$PI_NAME,
  TASK_ID = yaml$TASK_ID,
  PROJECT_LEAD_NAME = yaml$PROJECT_LEAD_NAME,
  DEPARTMENT = yaml$DEPARTMENT,
  LEAD_ANALYSTS = yaml$LEAD_ANALYSTS,
  GROUP_LEAD = yaml$GROUP_LEAD,
  CONTACT_EMAIL = yaml$CONTACT_EMAIL,
  PIPELINE = yaml$PIPELINE,
  START_DATE = yaml$START_DATE,
  COMPLETION_DATE = yaml$COMPLETION_DATE
)

#################################################################################
# Create reports_dir
results_dir <- file.path(analysis_dir, "reports")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)}

# Define dynamic output file name using Sys.Date()
output_file <- file.path(results_dir, glue::glue("Updates-{PROJECT_NAME}-{Sys.Date()}.html"))


# Render the R Markdown file
render("01-generate-project-report.Rmd", output_file = output_file, params = report_params)
#################################################################################

