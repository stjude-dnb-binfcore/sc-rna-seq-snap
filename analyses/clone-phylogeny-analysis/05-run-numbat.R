#################################################################################
# Load library
suppressPackageStartupMessages({
  library(knitr)
  library(tidyverse)
  library(yaml)
  library(optparse)
  library(Seurat)
  library(numbat)
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
analysis_dir <- file.path(root_dir, "analyses", "clone-phylogeny-analysis") 
module_results_dir <- file.path(analysis_dir, "results")
count_mat_results_dir <- file.path(analysis_dir, "results", "01-create-count-mat") 
df_allele_results_dir <- file.path(analysis_dir, "results", "04-create-df-allele-rda") 

# Input files
#reference_file <- file.path(input_dir, "myref-Tcellls.rda")
count_mat_file <- dir(path = count_mat_results_dir,  pattern = ("count_mat.rds"), full.names = TRUE, recursive = TRUE)
df_allele_file <- dir(path = df_allele_results_dir,  pattern = ("df_allele.rda"), full.names = TRUE, recursive = TRUE)

# Create step_results_dir
step_results_dir <- 
  file.path(module_results_dir, paste0("05-run-numbat"))
if (!dir.exists(step_results_dir)) {
  dir.create(step_results_dir)
}

#######################################################
# Read reference file
#reference_obj <- readRDS(reference_file)

#######################################################
# Create list 
count_mat_list <- list()
df_allele_list <- list()

#######################################################
# Read metadata file and define `sample_name`
metadata_file <- file.path(metadata_dir, "project_metadata.tsv") # metadata input file

# Read metadata file and define `sample_name`
project_metadata <- read.csv(metadata_file, sep = "\t", header = TRUE)
sample_name <- unique(as.character(project_metadata$ID))
sample_name <- sort(sample_name, decreasing = FALSE)
print(sample_name)

for (i in seq_along(sample_name)){
  
  # Create results_dir per sample
  results_dir <- file.path(step_results_dir, sample_name[i])
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)}

  # Read data per sample
  cat("Processing sample:", sample_name[i], "\n")

  # Read the count matrix
  cat("Reading the count matrix for sample:", sample_name[i], "\n")
  cat("Processing count_mat_list[[", i, "]]:", count_mat_file[i], "\n")
  count_mat_list[[i]] <- readRDS(count_mat_file[i]) # Use double brackets to assign to the list

  # Read the allele data
  cat("Reading the allele data for sample:", sample_name[i], "\n")
  cat("Processing df_allele_list[[", i, "]]:", df_allele_file[i], "\n")
  df_allele_list[[i]] <- readRDS(df_allele_file[i]) # Use double brackets to assign to the list
  
  # Rename the cells name to match the ones from the snap pipeline
  df_allele_list[[i]] <- df_allele_list[[i]] %>%
    mutate(cell = str_c(sample_name[i], ":", cell),
           cell = str_remove(cell, "-1$"))
  print(head(df_allele_list[[i]]))
  
  
  ###############################################################################################################
  # Debugging step
  #cat("Input check - count_mat:", is.null(count_mat_list[[i]]), "df_allele:", is.null(df_allele_list[[i]]), "\n")
  #Input check - count_mat: FALSE df_allele: FALSE 
  #cat("Processing sample:", sample_name[i], "Count matrix size:", dim(count_mat_list[[i]]), "Allele data size:", dim(df_allele_list[[i]]), "\n")
  
  # Print out the results
  #cat("Common cells:", length(common_cells), "\n")
  #cat("Missing in count_mat:", length(missing_in_count_mat), "\n")
  #cat("Missing in df_allele:", length(missing_in_df_allele), "\n")
  
  ###############################################################################################################  
  
  
  ############## ############## ############## ##############
  # Sanity check
  cat("Processing sample:", sample_name[i], "\n")
  cat("count_mat_list sample:", count_mat_file[i], "\n")
  cat("df_allele_list sample:", df_allele_file[i], "\n")
  
  cat("Length of sample_name:", length(sample_name), "\n")
  cat("Length of count_mat_list:", length(count_mat_list), "\n")
  cat("Length of df_allele_list:", length(df_allele_list), "\n")
  
  cat("Number of cells in count_mat_list[[i]]:", ncol(count_mat_list[[i]]), "\n")
  cat("Number of genes in count_mat_list[[i]]:", nrow(count_mat_list[[i]]), "\n")
  cat("Number of cells in df_allele_list[[i]]:", nrow(df_allele_list[[i]]), "\n") 
  cat("Number of genes in df_allele_list[[i]]:", length(unique(df_allele_list[[i]]$gene)), "\n")
  #print(head(df_allele_list[[i]]))
  ############## ############## ############## ##############
  
  # Run numbat
  cat("Running numbat for sample:", sample_name[i], "\n")
  out = run_numbat(
    count_mat = count_mat_list[[i]], # gene x cell integer UMI count matrix 
    lambdas_ref = ref_hca, #reference_obj; # reference expression profile, a gene x cell type normalized expression level matrix
    df_allele = df_allele_list[[i]], # allele dataframe generated by pileup_and_phase script
    #gtf = gtf_hg38, # provided upon loading the package
    #genetic_map = genetic_map_hg38, # provided upon loading the package
    genome = "hg38",
    t = 1e-04, # 1e-3
    tau = 0.7,
    max_iter = 2,
    #init_k = 3,
    ncores = 64, #10
    plot = TRUE,
    max_entropy = 0.8,
    min_cells = 100, #20
    min_LLR = 1, #5; No CNV remains after filtering by LLR in pseudobulks. Consider reducing min_LLR.
    skip_nj=TRUE,
    #segs_loh = segs_loh[[i]],
    #gamma = 5, # Overdispersion in allele counts (default: 20)
    #multi_allelic = TRUE, # by default this is FALSE
    #check_convergence = TRUE, # by default tzhis is FALSE
    out_dir = results_dir) # Directory for output results
  cat("Running numbat complete for sample:", sample_name[i], "\n")
}
