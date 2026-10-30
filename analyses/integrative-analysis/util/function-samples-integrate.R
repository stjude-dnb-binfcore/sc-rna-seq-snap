############################################################################################################
#' Function to integrate samples using Seurat
#' @param seurat_obj_list
#' @param nfeatures_value
#' @param num_dim_seurat
#' @param num_dim_seurat_integration
#' @param big_data
#' @param reference_list
#' @param Genome
#' @param Regress_Cell_Cycle
#' @param PCA_Feature_List
#'
#' @return
#' @export
#'
#' @examples
#' 
#'
seurat_integration <- function(seurat_obj_list, nfeatures_value, num_dim_seurat, num_dim_seurat_integration, 
                               big_data, reference_list, Genome, Regress_Cell_Cycle, PCA_Feature_List){
  
  set.seed(1234) # Make code reproducible
  
  seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x){
    
    cat("Conducting SCTransform for sample:", unique(x$ID), "\n")
    # the following command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
    # Transformed data will be available in the SCT assay, which is set as the default after running sctransform
    x <- SCTransform(x, variable.features.n = nfeatures_value, verbose = TRUE, return.only.var.genes = FALSE, conserve.memory = TRUE)}) #vst.flavor = "v2"
  
  cat("Selecting Intgration Features\n")
  features <- SelectIntegrationFeatures(object.list = seurat_obj_list, nfeatures = nfeatures_value, verbose = TRUE)
  seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = features, verbose = TRUE)
  seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x){
    x <- RunPCA(x, features = features, verbose = TRUE)})
  
  # We will consider projects that contain big data, i.e., >= 1.3M cells.
  if (big_data == TRUE){
    if (missing(reference_list)){
      cat("Finding Integration Anchors using rPCA reduction, but no Reference\n")
      anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT", anchor.features = features, dims = 1:num_dim_seurat_integration, reduction = "rpca", k.anchor = 20, verbose = TRUE)
      cat("Performing Seurat Integration\n")
      #seurat_obj <- IntegrateData(anchorset = anchors, new.assay.name = "SCT", normalization.method = "SCT", dims = 1:num_dim_seurat_integration, verbose = TRUE)
      seurat_obj <- IntegrateData(anchorset = anchors)
      
      } else {
        cat("Finding Integration Anchors using rPCA reduction and a Reference\n")
        anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, reference = reference_list, normalization.method = "SCT", anchor.features = features, dims = 1:num_dim_seurat_integration, reduction = "rpca", k.anchor = 20, verbose = TRUE)
        cat("Performing Seurat Integration\n")
        #seurat_obj <- IntegrateData(anchorset = anchors, new.assay.name = "SCT", normalization.method = "SCT", dims = 1:num_dim_seurat_integration, verbose = TRUE)
        seurat_obj <- IntegrateData(anchorset = anchors)}
    
    } else {
      if (missing(reference_list)) {
        cat("Finding Integration Anchors using default parameters\n")
        anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT", anchor.features = features, verbose = TRUE)
        #seurat_obj <- IntegrateData(anchorset = anchors, new.assay.name = "SCT", normalization.method = "SCT", verbose = TRUE)
        seurat_obj <- IntegrateData(anchorset = anchors)
        
        } else {
          cat("Finding Integration Anchors using a Reference\n")
          anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, reference = reference_list, normalization.method = "SCT", anchor.features = features, verbose = TRUE)
          #seurat_obj <- IntegrateData(anchorset = anchors, new.assay.name = "SCT", normalization.method = "SCT", verbose = TRUE)
          seurat_obj <- IntegrateData(anchorset = anchors)}
    }
 
  cat("Assigning Cell-Cycle Scores\n")
  seurat_obj <- Genome_Specific_Cell_Cycle_Scoring(seurat_obj, Genome)
  
  if (Regress_Cell_Cycle == "NO"){
    cat("Scaling the data\n")
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = TRUE)
    
    }else if (Regress_Cell_Cycle == "YES"){
      cat("Regressing out Cell Cycle\n")
      seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_obj), verbose = TRUE)
      
      } else if (Regress_Cell_Cycle == "DIFF"){
        cat("Regressing out the difference between the G2M and S phase scores\n")
        seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
        seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "CC.Difference", features = rownames(seurat_obj), verbose = TRUE)}
  
  if (missing(PCA_Feature_List)){
    cat("Performing linear dimensional reduction on highly variable genes\n")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = TRUE)
    } else{
      cat("Performing linear dimensional reduction on custom gene list\n")
      seurat_obj <- RunPCA(seurat_obj, features = PCA_Feature_List, verbose = TRUE)}
  
  cat("Performing Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique\n")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:num_dim_seurat, verbose = TRUE)
  
  cat("Generate and add metadata\n")
  reduction_names <- c(paste0("umap"), paste0("pca")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  seurat_obj@meta.data <- merge_metadata(seurat_obj, metadata) 
  write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_{integration_method}.tsv"))) # Save metadata
  
  cat("Computing the k.param nearest neighbors for a given dataset\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:num_dim_seurat, verbose = TRUE)
  
  cat("Calculating clusters\n")
  seurat_obj <- FindClusters(seurat_obj, verbose = TRUE)
  
  cat("Returned Seurat object")
  return(seurat_obj)}



############################################################################################################
#' Function to integrate samples using Harmony
#' @param seurat_obj
#' @param variables_to_integrate
#' @param num_dim
#'
#' @return
#' @export
#'
#' @examples
#' 
#'
harmony_integration <- function(seurat_obj, variables_to_integrate, num_dim){
  
  set.seed(1234) # Make code reproducible
  
  cat("Performing Harmony Integration\n")
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = variables_to_integrate, plot_convergence = TRUE)
  
  cat("Performing Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique\n")
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:num_dim)
  
  cat("Generate and add metadata\n")
  reduction_names <- c(paste0("umap"), paste0("pca")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  seurat_obj@meta.data <- merge_metadata(seurat_obj, metadata) 
  write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_{integration_method}.tsv"))) # Save metadata
  
  
  cat("Computing the k.param nearest neighbors for a given dataset and clusters\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:num_dim) %>%
    FindClusters() %>% 
    identity()
  
  cat("Returned Seurat object")
  return(seurat_obj)}
  

############################################################################################################
#' Function to integrate samples using Liger
#' @param seurat_obj
#' @param variables_to_integrate
#' @param n_neighbors
#'
#' @return
#' @export
#'
#' @examples
#' 
#'
liger_integration <- function(seurat_obj, variables_to_integrate, n_neighbors){
  
  set.seed(1234) # Make code reproducible
  
  seurat_obj <- seurat_obj %>%
    normalize() %>% # normalization by number of UMIs
    selectGenes() %>% # Variable gene selection
    scaleNotCenter() # scaling of individual genes but do not center gene expression because NMF requires nonnegative values
  
  cat("Performing Liger Integration\n")
  # Run integrative non-negative matrix factorization (iNMF) to learn a low-dimensional space 
  # in which each cell is defined by one set of dataset-specific factors, or metagenes, and another set of shared metagenes.
  # Each factor often corresponds to a biologically interpretable signal—like the genes that define a particular cell type. 
  # A tuning parameter, λ, allows adjusting the size of dataset-specific effects 
  # to reflect the divergence of the datasets being analyzed.
  # number of factors k
  
  # Identifying shared and dataset-specific factors through integrative non-negative matrix factorization (iNMF). 
  # We have derived and implemented a novel coordinate descent algorithm for efficiently performing the factorization.
  seurat_obj <- runIntegration(seurat_obj, k = 20) #lambda = 5, split.by = variables_to_integrate
  
  # After performing iNMF, we use a novel strategy that increases robustness of joint clustering. 
  # We first assign each cell a label based on the maximum factor loading, 
  # then build a shared factor neighborhood graph, in which we connect cells that have similar factor loading patterns.
  seurat_obj <- quantileNorm(seurat_obj) # `RunQuantileNorm` according to your needs; split.by = variables_to_integrate
  
  # Jointly clustering cells and normalizing factor loadings.
  # You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
  cat("Computing the k.param nearest neighbors\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "inmf", k.param = 20, dims = 1:n_neighbors)
  
  cat("Calculating clusters\n")
  seurat_obj <- FindClusters(seurat_obj)
  
  # Visualization using t-SNE or UMAP and analysis of shared and dataset-specific marker genes
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:ncol(seurat_obj[["inmf"]]), reduction = "inmf", verbose = TRUE)
  
  cat("Generate and add metadata\n")
  reduction_names <- c(paste0("umap"), paste0("pca")) # Export the reductions to Seurat
  metadata <- as_data_frame_seurat(seurat_obj, reduction = reduction_names, metadata = TRUE)
  seurat_obj@meta.data <- merge_metadata(seurat_obj, metadata) 
  write_tsv(metadata, file = paste0(results_dir, "/", glue::glue("metadata_{integration_method}.tsv"))) # Save metadata
  
  cat("Returned Seurat object")
  return(seurat_obj)}

  ############################################################################################################
  
  