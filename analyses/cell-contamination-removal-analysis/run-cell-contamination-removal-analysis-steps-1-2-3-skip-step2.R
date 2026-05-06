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
data_dir_module <- yaml$data_dir_module_name
umap_value <- yaml$umap_value

analysis_dir <- file.path(root_dir, "analyses", "cell-contamination-removal-analysis") 
report_dir <- file.path(analysis_dir, "plots") 

################################################################################################################
# Run Rmd scripts 
################################################################################################################
future_globals_value <- as.numeric(yaml$future_globals_value_contamination) * 1024^3
resolution = yaml$resolution_clustering_module
#integration_method = yaml$integration_method


##########################################################################################
# Define dual index genome names
dual_genomes_upstream <- c("DualGRCh38", "Dualhg19", "DualGRCm39", "Dualmm10", "Dualmm9")

# Select genome_name based on upstream genome type
if (yaml$genome_name_upstream %in% dual_genomes_upstream) {
  genome_name <- yaml$genome_name_contamination
  source_genome <- "contamination"
  cat("Using genome name: ", genome_name, "as defined in the contamination module\n")
  
} else {
  genome_name <- yaml$genome_name_upstream
  source_genome <- "upstream"
  cat("Using genome name: ", genome_name, "as defined in the upstream module\n")
}

##########################################################################################


# STEP 1 - Remove contamination
rmarkdown::render('01-cell-contamination-removal.Rmd', clean = TRUE,
                  output_dir = file.path(report_dir),
                  output_file = c(paste('Report-', 'cell-contamination-removal', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                # this will be the single resolution that fits the data the best
                                resolution_list = yaml$resolution_list_find_markers, 
                                # number of clusters to keep
                                keep_clusters = yaml$keep_clusters_contamination_module, 
                                
                                grouping = yaml$grouping,
                                min_genes = yaml$min_genes, 
                                min_count = yaml$min_count,
                                mtDNA_pct_default = yaml$mtDNA_pct_default,
                                
                                # process_seurat
                                nfeatures_value = yaml$nfeatures_value,
                                #genome_name = yaml$genome_name_upstream,
                                Regress_Cell_Cycle_value = yaml$Regress_Cell_Cycle_value,
                                assay = yaml$assay_contamination_module,
                                num_pcs = yaml$num_pcs,
                                prefix = yaml$prefix,
                                num_dim = yaml$num_dim_filter_object,
                                num_neighbors = yaml$num_neighbors_filter_object,
                                use_condition_split = yaml$use_condition_split_filter_object,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                print_pdf = yaml$print_pdf_filter_object,
                                PCA_Feature_List_value = yaml$PCA_Feature_List_value,
                                
                                # the following parameters are defined in the `yaml` file
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
# STEP 2 - Integration
# skip step 2 - single-sample cohort

################################################################################################################
# step 3 - Clustering
resolution = yaml$resolution_clustering_module

rmarkdown::render('03-cluster-cell-calling.Rmd', clean = TRUE,
                  output_dir = file.path(report_dir),
                  output_file = c(paste('Report-', glue::glue("cluster-cell-calling-{resolution}"), '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                num_dim = yaml$num_dim_clustering_module, 
                                reduction_value = yaml$reduction_value_clustering_module,
                                resolution_list = yaml$resolution_list_clustering_module, 
                                resolution_list_default = yaml$resolution_list_default_clustering_module,
                                algorithm_value = yaml$algorithm_value_clustering_module, 
                                assay = yaml$assay_contamination_module,
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

