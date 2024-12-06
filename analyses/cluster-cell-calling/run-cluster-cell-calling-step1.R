#################################################################################
# This will run all scripts in the module
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)})

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
analysis_dir <- file.path(root_dir, "analyses", "cluster-cell-calling") 
report_dir <- file.path(analysis_dir, "plots") 

################################################################################################################
# Run Rmd scripts 
################################################################################################################
future_globals_value = yaml$future_globals_value_clustering_module
resolution = yaml$resolution_clustering_module

rmarkdown::render('01-cluster-cell-calling.Rmd', clean = TRUE,
                  output_dir = file.path(report_dir),
                  output_file = c(paste('Report-', glue::glue("cluster-cell-calling-{resolution}"), '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                num_dim = yaml$num_dim_clustering_module, 
                                reduction_value = yaml$reduction_value_clustering_module,
                                resolution_list = yaml$resolution_list_clustering_module, 
                                resolution_list_default = yaml$resolution_list_default_clustering_module,
                                algorithm_value = yaml$algorithm_value_clustering_module, 
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
