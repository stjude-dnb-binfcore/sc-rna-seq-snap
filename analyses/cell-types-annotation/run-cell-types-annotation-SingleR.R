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
  library(celldex)
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

broad_report_dir <- file.path(module_plots_dir, "01_cell_types_annotation_SingleR_broad") 
if (!dir.exists(broad_report_dir)) {
  dir.create(broad_report_dir)}

fine_report_dir <- file.path(module_plots_dir, "02_cell_types_annotation_SingleR_fine") 
if (!dir.exists(fine_report_dir)) {
  dir.create(fine_report_dir)}

#gene_markers_report_dir <- file.path(module_plots_dir, "03_cell_types_annotation_gene_markers") 
#if (!dir.exists(fine_report_dir)) {
#  dir.create(fine_report_dir)}

################################################################################################################
# Celldex reference for SingleR (configured in project_parameters.Config.yaml as
# `singler_celldex_reference`: the celldex *function name* only, e.g. HumanPrimaryCellAtlasData).
# The reference object cannot live in YAML; we load it by name from the celldex package namespace.
# https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html
################################################################################################################
singler_celldex_allowlist <- c(
  "HumanPrimaryCellAtlasData",
  "BlueprintEncodeData",
  "DatabaseImmuneCellExpressionData",
  "ImmGenData",
  "MonacoImmuneData",
  "MouseRNAseqData",
  "NovershternHematopoieticData"
)
ref_name <- yaml$singler_celldex_reference
if (is.null(ref_name) || !nzchar(trimws(as.character(ref_name)))) {
  ref_name <- "HumanPrimaryCellAtlasData"
} else {
  ref_name <- trimws(as.character(ref_name))
}
if (!ref_name %in% singler_celldex_allowlist) {
  stop(
    "Invalid singler_celldex_reference: \"", ref_name, "\".\n",
    "Use one of: ", paste(singler_celldex_allowlist, collapse = ", "),
    "\nOr extend `singler_celldex_allowlist` in run-cell-types-annotation-SingleR.R."
  )
}
celldex_loader <- getExportedValue("celldex", ref_name)
if (!is.function(celldex_loader)) {
  stop("celldex::", ref_name, " is not a function; check your Bioconductor celldex version.")
}
bpe <- celldex_loader()

################################################################################################################

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
input_data_file <- file.path(data_dir_annotation_module, glue::glue("seurat_obj_clusters_all.rds"))

################################################################################################################
# Broad cell type annotation
rmarkdown::render('01-cell-types-annotation-SingleR-broad.Rmd', clean = TRUE,
                  output_dir = file.path(broad_report_dir),
                  output_file = c(paste('Report-', 'cell-types-annotation-SingleR-broad', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                reduction_value = yaml$reduction_value_annotation_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                min.diff.med_value = yaml$min.diff.med_value_annotation_module,
                                use_min.diff.med = yaml$use_min.diff.med_annotation_module,
                                data_file = input_data_file,
                                assay = yaml$assay_annotation_module,
                                ct_palette_file = yaml$ct_palette_file_broad,
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
# Fine cell type annotation
rmarkdown::render('02-cell-types-annotation-SingleR-fine.Rmd', clean = TRUE,
                  output_dir = file.path(fine_report_dir),
                  output_file = c(paste('Report-', 'cell-types-annotation-SingleR-fine', '-', Sys.Date(), sep = '')),
                  output_format = 'all',
                  params = list(integration_method = yaml$integration_method_clustering_module,
                                reduction_value = yaml$reduction_value_annotation_module,
                                condition_value1 = yaml$condition_value1,
                                condition_value2 = yaml$condition_value2,
                                condition_value3 = yaml$condition_value3,
                                min.diff.med_value = yaml$min.diff.med_value_annotation_module,
                                use_min.diff.med = yaml$use_min.diff.med_annotation_module,
                                data_file = input_data_file,
                                assay = yaml$assay_annotation_module,
                                ct_palette_file = yaml$ct_palette_file_fine,
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