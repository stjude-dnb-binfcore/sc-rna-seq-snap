############################################################################################################
#' Resolve de-go-analysis parameters when Rmd is knit without run-de-go-analysis.R
############################################################################################################

is_dego_param_placeholder <- function(x) {
  is.null(x) ||
    length(x) != 1L ||
    is.na(x) ||
    !nzchar(trimws(as.character(x))) ||
    identical(as.character(x), ".")
}

find_project_config_path <- function(root_dir = ".") {
  candidates <- unique(c(
    file.path(root_dir, "project_parameters.Config.yaml"),
    "../../project_parameters.Config.yaml",
    "../../../project_parameters.Config.yaml"
  ))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    return(NULL)
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

read_dego_project_yaml <- function(root_dir = ".") {
  cfg_path <- find_project_config_path(root_dir)
  if (is.null(cfg_path)) {
    return(NULL)
  }
  read_dego_config_yaml(cfg_path)
}

read_dego_config_yaml <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  lines <- readLines(config_path, encoding = "UTF-8", warn = FALSE)
  if (length(lines) == 0L) {
    stop("Configuration file is empty: ", config_path)
  }
  yaml::read_yaml(text = paste(lines, collapse = "\n"))
}

resolve_dego_string_param <- function(param_value,
                                      yaml_key,
                                      root_dir,
                                      label = yaml_key) {
  if (!is_dego_param_placeholder(param_value)) {
    return(as.character(param_value))
  }
  yaml_proj <- read_dego_project_yaml(root_dir)
  if (!is.null(yaml_proj) && !is.null(yaml_proj[[yaml_key]])) {
    val <- yaml_proj[[yaml_key]]
    if (!is_dego_param_placeholder(val)) {
      return(as.character(val))
    }
  }
  stop(
    label, " is not set. Use bash run-de-go-analysis.sh, or set ", yaml_key,
    " in project_parameters.Config.yaml (and pass root_dir in report params)."
  )
}

assert_seurat_metadata_column <- function(seurat_obj, column, context = "metadata") {
  if (!column %in% colnames(seurat_obj@meta.data)) {
    stop(
      "Column '", column, "' was not found in Seurat ", context, ". ",
      "Check cell_type_label_go_module in project_parameters.Config.yaml. ",
      "Available columns include: ",
      paste(utils::head(colnames(seurat_obj@meta.data), 40L), collapse = ", "),
      if (ncol(seurat_obj@meta.data) > 40L) ", ..." else ""
    )
  }
  invisible(column)
}

############################################################################################################
