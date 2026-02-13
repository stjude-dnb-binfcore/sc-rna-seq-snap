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
analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis") 

module_plots_dir <- file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)}

module_results_dir <- file.path(analysis_dir, "results") 
if (!dir.exists(module_results_dir)) {
  dir.create(module_results_dir)}
################################################################################################################

# Step 1
rmarkdown::render('01-de-analysis.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'de-analysis', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(assay = yaml$assay_annotation_module,
                                cell_type_label_go = yaml$cell_type_label_go_module,
                                resolution_list = yaml$resolution_list_find_markers, 
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                
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

################################################################################################################
# Step 2
rmarkdown::render('02-volcano-plots.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'volcano-plots', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(n_value = yaml$n_value_find_markers,
                                resolution_list = yaml$resolution_list_find_markers, 
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                
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

################################################################################################################
# Step 3
rmarkdown::render('03-go-clusterprofiler.Rmd', clean = TRUE,
                  output_dir = file.path(module_plots_dir),
                  output_file = c(paste('Report-', 'go-clusterprofiler', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(OrgDb_value = yaml$OrgDb_value,
                                resolution_list = yaml$resolution_list_find_markers, 
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                universe_value = yaml$universe_value_module,
                                show_number_pathways = yaml$show_number_pathways_value_module,

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



