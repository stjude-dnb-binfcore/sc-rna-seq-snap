#################################################################################
# This will run all scripts in the module
#################################################################################
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
})

#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)) {
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")
}

config_lines <- readLines(configFile, encoding = "UTF-8", warn = FALSE)
yaml <- read_yaml(text = paste(config_lines, collapse = "\n"))

if (is.null(yaml$cell_type_label_go_module) || !nzchar(as.character(yaml$cell_type_label_go_module))) {
  stop(
    "cell_type_label_go_module is missing or empty in ", configFile,
    ". Set it to a metadata column from cell-types-annotation (e.g. predicted.cell.signature.ident)."
  )
}

#################################################################################
# Set up directories and paths to root_dir and analysis_dir
root_dir <- yaml$root_dir
analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis")

module_plots_dir <- file.path(analysis_dir, "plots")
if (!dir.exists(module_plots_dir)) {
  dir.create(module_plots_dir)
}

module_results_dir <- file.path(analysis_dir, "results")
if (!dir.exists(module_results_dir)) {
  dir.create(module_results_dir)
}

source("util/function-de-global-interaction.R")
source("util/function-de-go-config.R")
dego_mode <- normalize_dego_comparison_mode(yaml$dego_condition_comparison_mode)
dego_plots_dir <- dego_module_plots_dir(analysis_dir, dego_mode)
if (!dir.exists(dego_plots_dir)) {
  dir.create(dego_plots_dir, recursive = TRUE)
}

future_globals_value <- as.numeric(yaml$future_globals_value_dego) * 1024^3

common_report_params <- list(
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
  COMPLETION_DATE = yaml$COMPLETION_DATE
)

if (run_dego_mode_needs_precompute(dego_mode)) {
  source("run-de-computation.R")
  precomputed <- run_de_go_computation(config_path = configFile, module_dir = ".")
  markers_all <- file.path(
    dego_results_de_dir(analysis_dir, dego_mode),
    paste0("Res_", yaml$resolution_list_find_markers, "_Markers_all.tsv")
  )
  if (!isTRUE(precomputed) || !file.exists(markers_all)) {
    stop(
      "DE pre-computation did not complete for dego_condition_comparison_mode=",
      dego_mode,
      ". Expected file missing: ", markers_all,
      ". Check job logs and project_parameters.Config.yaml encoding (UTF-8)."
    )
  }
}

################################################################################################################
# step 1
rmarkdown::render(
  "01-de-analysis.Rmd",
  clean = TRUE,
  output_dir = dego_plots_dir,
  output_file = paste("Report-", "de-analysis", "-", Sys.Date(), sep = ""),
  output_format = "all",
  params = c(
    list(
      assay = yaml$assay_annotation_module,
      cell_type_label_go = resolve_dego_string_param(
        yaml$cell_type_label_go_module,
        yaml_key = "cell_type_label_go_module",
        root_dir = yaml$root_dir,
        label = "cell_type_label_go_module"
      ),
      resolution_list = yaml$resolution_list_find_markers,
      dego_condition_comparison_mode = dego_mode,
      min_cells_per_group = yaml$min_cells_per_group_dego_module,
      dego_ident1 = yaml$dego_ident1,
      dego_ident2 = yaml$dego_ident2,
      future_globals_value = future_globals_value,
      future_globals_value_dego_gib = yaml$future_globals_value_dego,
      condition_value1 = yaml$condition_value1,
      condition_value2 = yaml$condition_value2,
      condition_value3 = yaml$condition_value3
    ),
    common_report_params
  )
)

################################################################################################################
# step 2
rmarkdown::render(
  "02-volcano-plots.Rmd",
  clean = TRUE,
  output_dir = dego_plots_dir,
  output_file = paste("Report-", "volcano-plots", "-", Sys.Date(), sep = ""),
  output_format = "all",
  params = c(
    list(
      n_value = yaml$n_value_find_markers,
      resolution_list = yaml$resolution_list_find_markers,
      dego_condition_comparison_mode = dego_mode,
      min_cells_per_group = yaml$min_cells_per_group_dego_module,
      condition_value1 = yaml$condition_value1,
      condition_value2 = yaml$condition_value2,
      condition_value3 = yaml$condition_value3
    ),
    common_report_params
  )
)

################################################################################################################
# step 3
rmarkdown::render(
  "03-go-clusterprofiler.Rmd",
  clean = TRUE,
  output_dir = dego_plots_dir,
  output_file = paste("Report-", "go-clusterprofiler", "-", Sys.Date(), sep = ""),
  output_format = "all",
  params = c(
    list(
      OrgDb_value = yaml$OrgDb_value,
      resolution_list = yaml$resolution_list_find_markers,
      dego_condition_comparison_mode = dego_mode,
      min_cells_per_group = yaml$min_cells_per_group_dego_module,
      condition_value1 = yaml$condition_value1,
      condition_value2 = yaml$condition_value2,
      condition_value3 = yaml$condition_value3,
      universe_value = yaml$universe_value_module,
      show_number_pathways = yaml$show_number_pathways_value_module
    ),
    common_report_params
  )
)
################################################################################################################