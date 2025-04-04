#################################################################################
# This will run all scripts in the module
#################################################################################
# Load the Package with a Specific Library Path
#.libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.4")
#################################################################################
# Load library
suppressPackageStartupMessages({
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
# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "cell-types-annotation") 

module_plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)}

gene_markers_report_dir <- file.path(module_plots_dir, "04_cell_types_annotation_reference") 
if (!dir.exists(gene_markers_report_dir)) {
  dir.create(gene_markers_report_dir)}

################################################################################################################

future_globals_value = 214748364800 # 200*1024^3; other options: 1000 * 1024^2 = 1048576000; 8000 * 1024^2 =8388608000
resolution = yaml$resolution_list_find_markers
resolution_for_input_data = yaml$resolution_find_markers
integration_method = yaml$integration_method_clustering_module
input_data = yaml$module_with_input_data
input_data_folder= yaml$input_data_folder_name

################################
# Set data_dir
# Caution! Sometimes this file will be located in the `cluster-cell-calling` module
# BUT if we had to remove contamination, then it will be located in the `cell-contamination-removal-analysis` module
data_dir_annotation_module <- file.path(root_dir, "analyses", input_data, "results", glue::glue("{input_data_folder}_{resolution_for_input_data}"))
input_data_file <- file.path(data_dir_annotation_module, glue::glue("seurat_obj_integrated_{integration_method}_clusters_all.rds"))

################################################################################################################
# Gene score cell type annotation
rmarkdown::render('04-cell-types-annotation-reference.Rmd', clean = TRUE,
                  output_dir = file.path(gene_markers_report_dir),
                  output_file = c(paste('Report-', 'cell-types-annotation-reference', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                redution_value = yaml$redution_value_annotation_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                min.diff.med_value = yaml$min.diff.med_value_annotation_module,
                                use_min.diff.med = yaml$use_min.diff.med_annotation_module,
                                data_file = input_data_file,
                                assay = yaml$assay_annotation_module,
                                ct_palette_file = yaml$ct_palette_file_reference,

                                reference_dir = yaml$reference_dir_annotation_module,
                                reference_file_name = yaml$reference_file_name_annotation_module,
                                genome_name = yaml$genome_name_upstream,
                                gene_name_convention_update = yaml$gene_name_convention_update_annotation_module,
                                resolution = yaml$resolution_list_find_markers,
                                celltype_reference = yaml$celltype_reference_module,

                                root_dir = yaml$root_dir,
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
                                COMPLETION_DATE = yaml$COMPLETION_DATE))

################################################################################################################