############################################################################################################
#' Run GO enrichment per cluster using closest genes
#'
#' @param closest_genes_all Dataframe of genes with a 'cluster' and 'gene_id' column
#' @param OrgDb_value Annotation database (e.g., org.Hs.eg.db)
#' @param plots_dir Directory to save plots
#' @param pval_cutoff Adjusted p-value cutoff (default = 0.05)
#' @param qval_cutoff q-value cutoff (default = 0.05)
#' @param ontology GO ontology: "BP", "MF", or "CC"
#' @param key_type Gene ID type ("ENSEMBL", "SYMBOL", etc.)
#' @param universe_value 
#' @param show_number_pathways
#'
#' @return A named list of enrichGO results (class: enrichResult)
#' @export
#'
#' @examples
#' 
run_clusterwise_GO_enrichment_markers <- function(closest_genes_all,
                                          OrgDb_value,
                                          plots_dir,
                                          pval_cutoff = 0.05,
                                          qval_cutoff = 0.05,
                                          ontology = "BP",
                                          key_type = "ENSEMBL",
                                          universe_value = universe_value,
                                          show_number_pathways = show_number_pathways) {

  ego_results <- list()
  
  for (clust in sort(unique(closest_genes_all$cluster))) {
    message("Running GO enrichment for: ", clust)
    
    genes <- closest_genes_all %>%
      filter(cluster == clust) %>%
      pull(gene) %>%
      unique()
    
    if (length(genes) == 0) {
      warning("No gene IDs for ", clust)
      next
    }
    
    # Run enrichGO
    ego <- tryCatch({
      enrichGO(
        gene          = genes,
        universe      = universe_value, # names(geneList)
        #keyType       = key_type,
        keyType = "SYMBOL",
        OrgDb         = OrgDb_value,
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff,
        qvalueCutoff  = qval_cutoff,
        readable      = TRUE
      )
    }, error = function(e) {
      warning("enrichGO failed for ", clust, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(ego) || nrow(ego) == 0) {
      message("No significant GO terms for: ", clust)
      next
    }
    
    # Add pairwise term similarity (for emapplot)
    ego <- pairwise_termsim(ego)
    ego_results[[as.character(clust)]] <- ego
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      # Barplot
      cat("Attempting barplot for:", clust, "\n")
      try({
        suppressWarnings({
          bar <- print(barplot(ego, showCategory = show_number_pathways) +
            ggtitle(glue::glue("GO Enrichment - {clust} (Barplot)")))
          ggsave(file.path(plots_dir, glue::glue("barplot-{clust}-ego.png")),
                 plot = bar, width = 15, height = 8)
        })
      }, silent = TRUE)
      
      # Dotplot
      cat("Attempting dotplot for:", clust, "\n")
      try({
        suppressWarnings({
          dot <- print(dotplot(ego, showCategory = show_number_pathways) +
            ggtitle(glue::glue("GO Enrichment - {clust} (Dotplot)")))
          ggsave(file.path(plots_dir, glue::glue("dotplot-{clust}-ego.png")),
                 plot = dot, width = 15, height = 8)
        })
      }, silent = TRUE)
      
      # Emapplot
      cat("Attempting emapplot for:", clust, "\n")
      try({
        suppressWarnings({
          emap <- print(emapplot(ego, showCategory = show_number_pathways, layout = "kk") +
            ggtitle(glue::glue("GO Enrichment - {clust} (Emapplot)")))
          ggsave(file.path(plots_dir, glue::glue("emapplot-{clust}-ego.png")),
                 plot = emap, width = 15, height = 8)
      })
    }, silent = TRUE)
      
      } else {
        message("No significant GO terms for: ", clust)
        }
    
  }
  message("Saved GO enrichment plots for ", clust)
  return(ego_results)
}
############################################################################################################
