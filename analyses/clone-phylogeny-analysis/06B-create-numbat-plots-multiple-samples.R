#################################################################################
# This will run all scripts in the module
#################################################################################
# Load library
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(Seurat)
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
metadata_dir <- yaml$metadata_dir
analysis_dir <- file.path(root_dir, "analyses", "clone-phylogeny-analysis") 
module_results_dir <- file.path(analysis_dir, "results")
annotation_results_dir <- file.path(root_dir, "analyses", "cell-types-annotation", "results") 
annotations_dir <- yaml$annotations_dir_rshiny_app
annotations_all_results_dir <- file.path(annotation_results_dir, annotations_dir) 

# Input files
seurat_obj_file <- dir(path = annotations_all_results_dir, pattern =  "seurat_obj.*\\.rds", full.names = TRUE, recursive = TRUE)
print(seurat_obj_file)


# File path to plots directory 
module_plots_dir <-
  file.path(analysis_dir, "plots") 
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)
}

numbat_plots_dir <-
  file.path(module_plots_dir, "06-create-numbat-plots") 
if (!dir.exists(numbat_plots_dir)) {
  dir.create(numbat_plots_dir)
}


# Create numbat_dir
# numbat_results_dir <- 
#  file.path(module_results_dir, paste0("06-create-numbat-plots"))
#if (!dir.exists(numbat_results_dir)) {
#  dir.create(numbat_results_dir)
#}


#######################################################
# Read metadata file and define `sample_name`
#df_allele_results_dir <- file.path(module_results_dir, "03-create-df-allele-rda") 
#df_allele_file <- c(dir(path = df_allele_results_dir,  pattern = "df_allele.rda", full.names = TRUE, recursive = TRUE))

################################################################################################################
# Create plots
#sample_name <- c()

#for (i in seq_along(df_allele_file)) {
  # Extract sample name from the file path
#  sample_name <- c(sample_name, gsub("Create-", "", str_split_fixed(df_allele_file[i], "/", 16)[, 15]))
#}

# Sort the sample names outside the loop if needed
#sample_name <- sort(sample_name, decreasing = FALSE)
#print(sample_name)  # Print all sample names to check if they are stored correctly


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
# You do not need to repeat normalization or UMAP after splitting the Seurat object by sample. 
# The normalization will carry over, and the UMAP embeddings will be available in each individual Seurat object that you split.
seurat_obj_list <- SplitObject(seurat_obj, split.by = "orig.ident")
#seurat_obj_list

# Check if the number of samples matches the number of Seurat objects
if (length(seurat_obj_list) != length(sample_name)) {
  stop("Number of samples in the Seurat object list does not match the number of unique sample names.")
}

##########################################################################################################

for (i in seq_along(sample_name)){
  
  # Read data per sample
  #cat("Processing for sample:", sample_name[i], "\n")

  # seurat_obj
  #seurat_obj <- seurat_obj_list[[i]]

  # Create directory to save html reports
  samples_plots_dir <- file.path(numbat_plots_dir, sample_name[i]) 
  if (!dir.exists(samples_plots_dir)) {
    dir.create(samples_plots_dir)}
  
  # Create results_dir per sample
  #results_dir <- file.path(numbat_results_dir, sample_name[i])
  #if (!dir.exists(results_dir)) {
  #  dir.create(results_dir)}
  
  
  rmarkdown::render('06A-create-numbat-plots.Rmd', clean = TRUE,
                  output_dir = file.path(samples_plots_dir),
                  output_file = c(paste('Report-', 'create-numbat-plots-', sample_name[i], '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(cell_type_label = yaml$cell_type_label_numbat,
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
  
}
################################################################################################################
