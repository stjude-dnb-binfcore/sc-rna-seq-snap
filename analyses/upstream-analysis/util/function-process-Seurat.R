#############################################################################################
#' Function to calculate cell cycle scores (S and G2M) for each cell.
#' Genome specifies which gene set to use for cell cycle genes
#' @param seurat_obj
#' @param Genome
#'
#' @return
#' @export
#'
#' @examples
#' 
Genome_Specific_Cell_Cycle_Scoring <- function(seurat_obj, Genome) {
  
  # Pull standard cell cycle genes for cell cycle scoring
  human.s.genes <- cc.genes.updated.2019$s.genes
  human.g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  # Had to append "hg19-" to the beginning of all genes due to dual index reference genome used
  dual.hg19.s.genes <- paste("hg19-", human.s.genes, sep = "")
  dual.hg19.g2m.genes <- paste("hg19-", human.g2m.genes, sep = "")
  
  # Had to append "GRCh38-" to the beginning of all genes due to dual index reference genome used
  dual.GRCh38.s.genes <- paste("GRCh38-", human.s.genes, sep = "")
  dual.GRCh38.g2m.genes <- paste("GRCh38-", human.g2m.genes, sep = "")
  
  # Changed genes to title case for mice genes (Human genes are annotated in all CAPS)
  mm10.s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
  mm10.g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
  
  # Appended "GRCh38-" to the beginning of all genes due to dual index reference genome used
  dual.mm10.s.genes <- paste("mm10---", mm10.s.genes, sep = "")
  dual.mm10.g2m.genes <- paste("mm10---", mm10.g2m.genes, sep = "")
  
  if (Genome == "GRCh38" | Genome == "hg19") {
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = human.s.genes, g2m.features = human.g2m.genes)
  } else if (Genome == "mm10") {
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = mm10.s.genes, g2m.features = mm10.g2m.genes)
  } else if (Genome == "Dualhg19") {
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = dual.hg19.s.genes, g2m.features = dual.hg19.g2m.genes)
  } else if (Genome == "DualGRCh38") {
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = dual.GRCh38.s.genes, g2m.features = dual.GRCh38.g2m.genes)
  } else if (Genome == "Dualmm10") {
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = dual.mm10.s.genes, g2m.features = dual.mm10.g2m.genes)
  }
}



#############################################################################################
#' Function to process and cluster cell for an individual sample
#' @param seurat_obj
#' @param nfeatures_value
#' @param Genome
#' @param Regress_Cell_Cycle 
#' @param assay
#' @param num_pcs
#' @param prefix
#' @param num_dim
#' @param num_neighbors
#' @param results_dir
#' @param plots_output
#' @param use_condition_split
#' @param condition
#' @param print_pdf
#' @param PCA_Feature_List
#'
#' @return
#' @export
#'
#' @examples
#' 
Process_Seurat <- function(seurat_obj, nfeatures_value, Genome, Regress_Cell_Cycle, assay, 
                          num_pcs, prefix, num_dim, num_neighbors, results_dir, 
                          plots_output, use_condition_split, condition, print_pdf, PCA_Feature_List) {
  
  set.seed(1234) # Make code reproducible

  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures_value) 
  
  cat("Get cell cycle scores and then scale data", "\n")
  # Get cell cycle scores (S and G2M) for each cell, Genome specifies which gene set to use for cell cycle genes
  seurat_obj <- Genome_Specific_Cell_Cycle_Scoring(seurat_obj, Genome)
  
  cat("Scale data", "\n")
  # Indicates whether or not to regress for cell cycle and, if so, which method to use and scale data
  if (Regress_Cell_Cycle == "NO") {
    print("We will NOT regress out for cell cycle and scale data!")
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
    } else if (Regress_Cell_Cycle == "YES") {    
      print("We will regress out for cell cycle and scale data!")
      seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_obj))
      } else if (Regress_Cell_Cycle == "DIFF") {
        print("We will regress out for cell cycle difference and scale data!")
        seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
        seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "CC.Difference", features = rownames(seurat_obj))
  } 
  
  cat("Run PCA", "\n") # Get PCs for dimensionality reduction
  # Either use the prespecified PCA.Feature.List or use Variable Features found before
  if (missing(PCA_Feature_List)) {
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), assay = assay, npcs = num_pcs, verbose = FALSE)
  } else {
      seurat_obj <- RunPCA(seurat_obj, features = PCA_Feature_List, assay = assay, npcs = num_pcs, verbose = FALSE)
  }
  
  cat("Run UMAP", "\n")
  set.seed(1234) # Make code reproducible
  
  # Create a dataframe of all of the possible combinations of number of PCs to use for UMAP, and the number of neighbors
  num_dim_vect <- c(num_dim)
  num_neighbors_vect <- c(num_neighbors)
  possibilities <- expand.grid(num_dim_vect, num_neighbors_vect)
  
  # For each of these combinations, let's calculate UMAP
  for(k in 1:nrow(possibilities)) {
    num_dim <- possibilities[k, 1]
    num_neighbors <- possibilities[k, 2]
    seurat_obj <- run_dr(data = seurat_obj, dr_method = "umap", reduction = paste0("pca"),
                         num_dim_use = num_dim, assay = "RNA", num_neighbors = num_neighbors,
                         prefix = glue::glue("ndim{num_dim}nn{num_neighbors}"))}
  
  # Generate metadata
  reduction_names <- c(paste0("umap", "ndim", possibilities[,1], "nn", possibilities[,2]), paste0("pca")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  
  # Save output files
  write_tsv(metadata, file = paste0(results_dir, "/", "metadata", ".tsv")) # Save metadata
  
  cat("Plot the plots", "\n")
  # Identify the 10 most highly variable genes
  print(top10 <- head(VariableFeatures(seurat_obj), 10))
  
  # plot variable features with and without labels
  #plot1 <- VariableFeaturePlot(seurat_obj)
  #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #print(plot1 + plot2)
  
  # PCA ########################################################################
  # Examine and visualize PCA results a few different ways
  print(seurat_obj[[paste0("pca")]], dims = 1:5, nfeatures = 5)
  print(VizDimLoadings(seurat_obj, dims = 1:2, reduction = paste0("pca")))
  
  # Heatmap
  # It can be useful when trying to decide which PCs to include for further downstream analyses. 
  # Both cells and features are ordered according to their PCA scores. 
  # Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, 
  # which dramatically speeds plotting for large datasets.
  print(DimHeatmap(seurat_obj, reduction = paste0("pca"), dims = 1, cells = 500, balanced = TRUE))
  print(DimHeatmap(seurat_obj, reduction = paste0("pca"), dims = 1:15, cells = 500, balanced = TRUE))
  
  # Determine the ‘dimensionality’ of the dataset
  print(ElbowPlot(seurat_obj, reduction = paste0("pca")))
  
  cat("Plot UMAP", "\n")
  umap_values <- c(paste0("UMAP", "ndim", possibilities[,1], "nn", possibilities[,2]))
  
  if (print_pdf == "YES") {
    for (umap_val in umap_values) {
      print(umap_val)
      set.seed(1234) # Make code reproducible for UMAPs
    
      # Save plots
      umap_out <- file.path(paste0(plots_output, "/", umap_val, "/"))
      dir.create(umap_out)
      print(umap_out)
    
      # nFeature_RNA
      name <- paste0(umap_out, "01_", umap_val, "_nFeature_RNA.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$nFeature_RNA,
                                      palette = gradient_palette_df,
                                      title_name = "nFeature_RNA")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # nCount_RNA
      name <- paste0(umap_out, "02_", umap_val, "_nCount_RNA.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$nCount_RNA,
                                      palette = gradient_palette_df,
                                      title_name = "nCount_RNA")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # percent.mito
      name <- paste0(umap_out, "03_", umap_val, "_percent.mito.pdf")
      p <- create_UMAP_color_gradient(df = metadata,
                                      umap_val = umap_val,
                                      color_value = metadata$percent.mito,
                                      palette = gradient_palette_df,
                                      title_name = "percent.mito")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # orig.ident
      name <- paste0(umap_out, "04_", umap_val, "_orig.ident.pdf")
      p <- create_UMAP_orig_ident(df = metadata,
                                  umap_val = umap_val,
                                  color_value = "orig.ident",
                                  title_name = "orig.ident")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
      # ID
      name <- paste0(umap_out, "05_", umap_val, "_ID.pdf")
      p <- create_UMAP_ID(df = metadata,
                          umap_val = umap_val,
                          color_value = "ID",
                          title_name = "ID")
      pdf(file = name, width = 6, height = 5)
      print(p)
      dev.off()
    
     # condition
     name <- paste0(umap_out, "06_", umap_val, "_condition.pdf")
     p <- create_UMAP_condition(df = metadata,
                                umap_val = umap_val,
                                color_value = condition)
     pdf(file = name, width = 6, height = 5)
     print(p)
     dev.off()
     
     if (use_condition_split == "YES") {
       print("Print condition_split!")
       # condition_split
       name <- paste0(umap_out, "07_", umap_val, "_condition_split.pdf")
       p <- create_UMAP_condition_split(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata[[condition]])
       pdf(file = name, width = 10, height = 6)
       print(p)
       dev.off()
       } else { 
         print("Do NOT Print condition_split!")
         next }
     }
    } else if (print_pdf == "NO") {
      for (umap_val in umap_values) {
        print(umap_val)
        set.seed(1234) # Make code reproducible for UMAPs
    
        # Save plots
        umap_out <- file.path(paste0(plots_output, "/", umap_val, "/"))
        dir.create(umap_out)
        print(umap_out)
    
        # nFeature_RNA
        name <- paste0(umap_out, "01_", umap_val, "_nFeature_RNA.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$nFeature_RNA,
                                        palette = gradient_palette_df,
                                        title_name = "nFeature_RNA")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # nCount_RNA
        name <- paste0(umap_out, "02_", umap_val, "_nCount_RNA.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$nCount_RNA,
                                        palette = gradient_palette_df,
                                        title_name = "nCount_RNA")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # percent.mito
        name <- paste0(umap_out, "03_", umap_val, "_percent.mito.png")
        p <- create_UMAP_color_gradient(df = metadata,
                                        umap_val = umap_val,
                                        color_value = metadata$percent.mito,
                                        palette = gradient_palette_df,
                                        title_name = "percent.mito")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # orig.ident
        name <- paste0(umap_out, "04_", umap_val, "_orig.ident.png")
        p <- create_UMAP_orig_ident(df = metadata,
                                    umap_val = umap_val,
                                    color_value = "orig.ident",
                                    title_name = "orig.ident")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # ID
        name <- paste0(umap_out, "05_", umap_val, "_ID.png")
        p <- create_UMAP_ID(df = metadata,
                            umap_val = umap_val,
                            color_value = "ID",
                            title_name = "ID")
        ggsave(file = name, width = 6, height = 5, device = "png")

        # condition
        name <- paste0(umap_out, "06_", umap_val, "_condition.png")
        p <- create_UMAP_condition(df = metadata,
                                   umap_val = umap_val,
                                   color_value = condition)
        ggsave(file = name, width = 6, height = 5, device = "png")

        if (use_condition_split == "YES"){
          print("Print condition_split!")
          # condition_split
          name <- paste0(umap_out, "07_", umap_val, "_condition_split.png")
          p <- create_UMAP_condition_split(df = metadata,
                                           umap_val = umap_val,
                                           color_value = condition)
          ggsave(file = name, width = 10, height = 6, device = "png")
          } else { 
            print("Do NOT Print condition_split!")
            next }
    }
  }
  
  cat("Returned Seurat object", "\n")
  return(seurat_obj)
}


