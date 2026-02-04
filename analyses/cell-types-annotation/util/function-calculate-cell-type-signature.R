#' Function to calculate cell type signature
#' @param seurat_obj
#' @param results_dir
#' @param plots_dir
#' @param gene_markers_df
#' @param genome_name
#' @param gene_name_convention_update
#' @param assay
#' @param resolution
#' 
#' @return
#' @export
#'
#' @examples
#' 
calculate_cell_type_signature <- function(seurat_obj, results_dir, plots_dir, gene_markers_df, genome_name, gene_name_convention_update,
                                          assay, resolution) {

  pdf(file = paste(plots_dir, "/", "all_cell_type_signatures.pdf", sep = ""), width = length(gene_markers_df)*6, height = 6)
  
  merged_gene.markers <- NULL  # Initialize merged_gene.markers outside the loop

  for (i in 1:ncol(gene_markers_df)) {

    # Create cell type name for the current column
    cell.type <- paste(colnames(gene_markers_df)[i], ".score", sep = "")
    
    # Extract non-empty gene markers for the current column
    gene.markers <- as.character(gene_markers_df[, i][gene_markers_df[, i] != ""])
    
      if (genome_name == "GRCh38" | genome_name == "hg19"){
        gene.markers <- toupper(gene.markers)
        
        } else if (genome_name == "GRCm39" | genome_name == "mm10" | genome_name == "mm9"){
          gene.markers <- str_to_title(gene.markers)

             } else if (genome_name == "Dualhg19"){
               gene.markers <- paste("hg19-", toupper(gene.markers), sep = "")
            
               } else if (genome_name == "DualGRCh38"){
                gene.markers <- paste("GRCh38-", toupper(gene.markers), sep = "")
              
                } else if (genome_name == "Dualmm10"){
                 gene.markers <- paste("mm10---", str_to_title(gene.markers), sep = "") 
                 
                } else if (genome_name == "Dualmm9"){
                  gene.markers <- paste("mm9---", str_to_title(gene.markers), sep = "") 
      
                  } else if (genome_name == "DualGRCm39") {
                   gene.markers <- paste("GRCm39---", str_to_title(gene.markers), sep = "") } 

        # Filter gene markers that are present in rownames of 'seurat_obj'
        gene.markers <- gene.markers[gene.markers %in% rownames(seurat_obj)] 
        
        # Create a data frame for gene markers with the correct column name
        gene.markers_df <- data.frame(gene = gene.markers)
        
        # Rename the first (and only) column to the cell type name
        colnames(gene.markers_df)[1] <- cell.type

        
        # If merged_gene.markers is not NULL, adjust row count by padding with NA
        if (!is.null(merged_gene.markers)) {
          # Get the maximum row count (length)
          max_rows <- max(nrow(merged_gene.markers), nrow(gene.markers_df))
          
          # Pad with NAs if necessary
          merged_gene.markers <- merge(merged_gene.markers, gene.markers_df, by = "row.names", all = TRUE) %>%
            select(-Row.names)

        } else {
          # For the first iteration, just assign the current gene.markers_df
          merged_gene.markers <- gene.markers_df 
        }
        
        # Save gene.markers used for the annotation
        write_tsv(gene.markers_df, file = paste0(results_dir, "/", cell.type, "-gene.markers", ".tsv")) 
        write_tsv(merged_gene.markers, file = paste0(results_dir, "/", "merged_gene.markers", ".tsv")) 
        seurat_obj <- AddModuleScore(seurat_obj, features = list(gene.markers), name = cell.type)
    
        # Print plots  
        print(VlnPlot(seurat_obj, features = paste(cell.type, "1", sep = ""), pt.size = 0))
        print(VlnPlot(seurat_obj, features = paste(cell.type, "1", sep = ""), pt.size = 0, group.by = "ID"))
        print(VlnPlot(seurat_obj, features = paste(cell.type, "1", sep = ""), pt.size = 0, group.by = glue::glue("{assay}_snn_res.{resolution}")))}
    
    cell_cluster_scores <- seurat_obj@meta.data[, grepl(".score1", colnames(seurat_obj@meta.data))]
  
    # Predict.Cell.Score <- function(Data.Frame)
    # {
    #   Data.Frame <- sort(Data.Frame, decreasing = TRUE)
    #   if (foldchange(Data.Frame[1], Data.Frame[2]) > 1) {
    #     celltype <- sub(".score1", "", names(Data.Frame[1]))
    #   } else {
    #     celltype <- "UNDETERMINED" }
    #   return(celltype)
    # }
    # seurat_obj$predicted.cell.signature.ident <- apply(cell_cluster_scores, MARGIN = 1, FUN = Predict.Cell.Score)
    
    seurat_obj$predicted.cell.signature.ident <- sub(".score1", "", colnames(cell_cluster_scores)[max.col(cell_cluster_scores, ties.method = "first")])
  
    # Print plots
    print(DimPlot(seurat_obj, reduction = reduction_value, group.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
    print(DimPlot(seurat_obj, reduction = reduction_value, group.by = "Phase") + theme(aspect.ratio = 1))
    print(DimPlot(seurat_obj, reduction = reduction_value, group.by = "ID") + theme(aspect.ratio = 1))
    print(DimPlot(seurat_obj, reduction = reduction_value, split.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
    print(ggplot(data = data.frame(table(seurat_obj$predicted.cell.signature.ident)), aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + geom_text(aes(label = Freq)))
    
    dev.off()
  
    # Save tables
    #write.table(table(seurat_obj$predicted.cell.signature.ident, seurat_obj$ID), file = paste(results_dir, "/", "all_cell_type_signatures_count_data.tsv", sep = ""), sep="\t", quote = FALSE)
    #write.table(prop.table(table(seurat_obj$predicted.cell.signature.ident, seurat_obj$ID), margin = 2) *100, file = paste(results_dir, "/", "all_cell_type_signatures_percentage_data.tsv", sep = ""), sep="\t", quote=FALSE)
  
    return(seurat_obj)
}

