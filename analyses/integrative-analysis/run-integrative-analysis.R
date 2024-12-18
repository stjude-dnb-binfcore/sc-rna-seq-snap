#################################################################################
# This will run all scripts in the module
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)})

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
analysis_dir <- file.path(root_dir, "analyses", "integrative-analysis") 
report_dir <- file.path(analysis_dir, "plots") 

################################################################################################################
# Run Rmd scripts to process data per method
################################################################################################################
integration_method=yaml$integration_method

rmarkdown::render('01-integrative-analysis.Rmd', clean = TRUE,
                  output_dir = file.path(report_dir),
                  output_file = c(paste('Report-', 'integrative-analysis-{integration_method}', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(
                    # the following parameters are defined in the `yaml` file
                    future_globals_value = yaml$future_globals_value_integration_module,
                    use_seurat_integration = yaml$use_seurat_integration,
                    use_harmony_integration = yaml$use_harmony_integration,
                    use_liger_integration = yaml$use_liger_integration,
                    integration_method = yaml$integration_method,
                    num_dim_seurat =yaml$num_dim_seurat,
                    num_dim_seurat_integration = yaml$num_dim_seurat_integration,
                    big_data_value = yaml$big_data_value, 
                    num_dim_harmony = yaml$num_dim_harmony,
                    n_neighbors_value = yaml$n_neighbors_value,
                    variable_value = yaml$variable_value,
                    reference_list_value = yaml$reference_list_value, 
                    PCA_Feature_List_value = yaml$PCA_Feature_List_value,       
                    genome_name = yaml$genome_name,
                    nfeatures_value = yaml$nfeatures_value,
                    Regress_Cell_Cycle_value = yaml$Regress_Cell_Cycle_value,
                    root_dir = yaml$root_dir,
                    metadata_dir = yaml$metadata_dir,
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
