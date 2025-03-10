# Create color palettes for `sc-rna-seq-snap` repo
# Antonia Chroni <antonia.chroni@stjude.org> for DNB Bioinformatics Core Analysis Team
#
# Usage:
# Anywhere a plot is being made, source these TSV file and use the color palette for
# each appropriate data type.
#
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(celldex)
  #library(RColorBrewer)
  library(viridis)
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

# Output to palette directory
output_dir <-
  file.path(root_dir, "figures", "palettes")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

########################################################################################################################
# Get the number of unique cell types
bpe <- celldex::HumanPrimaryCellAtlasData()

unique_cell_types <- unique(bpe$label.fine)
num_unique_cell_types <- length(unique_cell_types)

# Generate a viridis color palette for the number of unique cell types (e.g., 157)
palette <- viridis(num_unique_cell_types)

# Print the color palette
#print(palette)

# Create a data frame with cell types and their corresponding colors
cell_type_colors <- data.frame(
  cell_type_names = unique_cell_types,
  hex_codes = palette
)

# Save the data frame as a .tsv file
write_tsv(cell_type_colors, file = paste0(output_dir, "/", "HumanPrimaryCellAtlasData_cell_types_with_colors.tsv")) # Save df

# Print a message to confirm
cat("The cell types and their hex colors have been saved to 'unique_cell_types_with_colors.tsv'.\n")

