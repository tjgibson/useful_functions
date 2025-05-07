# Generate table of gene overlaps with one column per sample, rows as genes, and TRUE/FALSE for if each gene is present in a given sample
gene_overlap <- function(x, gene_col, sample_col) {
  out_table <- tibble(
    gene_id = sort(unique(x[[gene_col]]))
  )
  
  for (i in unique(x[[sample_col]])) {
    i_data <- x[(x[[sample_col]] == i),]
    out_table[[i]] <- out_table$gene_id %in% i_data[[gene_col]]
  }
  
  return(out_table)
}

# pairwise gene overlap matrix, returning either n genes or percent
# useful for making heatmaps of gene list similarity
pairwise_gene_overlap <- function(x, gene_col, sample_col, return_type = "n") {
  groups <- unique(x[[sample_col]])
  out_table <- matrix(nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
  
  if (!(return_type %in% c("n", "percent"))) {
    stop("return_type should be `n` or `percent`")
  }
  
  for (i in 1:length(groups)) {
    for (j in 1:length(groups)) {
      i_genes <- x[(x[[sample_col]] == groups[i]),][[gene_col]]
      j_genes <- x[(x[[sample_col]] == groups[j]),][[gene_col]]
      
      if (return_type == "n") {
        out_table[i,j] <- sum(i_genes %in% j_genes)
      }
      else if (return_type == "percent") {
        out_table[i,j] <- sum(i_genes %in% j_genes) / length(i_genes)
      }
    }
    
  }
  
  return(out_table)
}