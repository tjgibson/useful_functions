#functions for processing peaks (e.g. ChIP-seq, ATAC-seq) ----------------------
# function to extend peak summits of Granges object
extend_summits <- function(x, upstream=100L, downstream=100L) {
  require(GenomicRanges)
  
  # check input
  if (!(inherits(x, "GRanges") | inherits(x, "GRangesList"))) {
    stop("x must be a GRanges or GRangesList")
  }
  
  
  
  if (!is.integer(upstream)) {
    stop("upstream must be an integer")
  }
  
  if (!is.integer(downstream)) {
    stop("downstream must be an integer")
  }
  
  
  
  start(x) <- start(x) - abs(upstream)
  end(x) <- end(x) + abs(downstream)
  return(x)
}


# function to merge multiple sets of ChIP peaks into a non-overlapping master set
# supports multiple methods
# peaks should be a GRangesList
# method should be either "overlap" or "merge": 
# overlap method will determine overlaps among all peak sets and combine non-overlapping peaks
# merge method will combine overlapping peaks into a single peak
combine_peaks <- function(peaks, method = "overlap", min_overlap = 1L) {
  require(GenomicRanges)
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # merge peaks using "overlap" method
  if (method == "overlap") {
    combined_peaks <- peaks[[1]]
    for (i in 2:length(peaks)) {
      i_set <- peaks[[i]]
      combined_peaks <- c(combined_peaks,
                          subsetByOverlaps(i_set, combined_peaks, minoverlap = min_overlap, invert = TRUE))
    }
    
  }
  # merge peaks using "merge" method
  if (method == "merge") {
    combined_peaks <- GenomicRanges::reduce(unlist(peaks))
    
  }
  mcols(combined_peaks) <- NULL
  return(combined_peaks)
  
}

# function to build an overlap table
# takes a Granges list object
# combines all nonoverlapping peaks from all samples, then builds a logical table indicating which peaks from each sample overlap each peak on the master list
peak_overlap_table <- function(peaks, method = "overlap", min_overlap = 1L, combine_peaks = TRUE) {
  
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  
  # set names
  if(is.null(names(peaks))) {
    names(peaks) <- paste0("sample_", seq((peaks)))
  }
  
  
  if (combine_peaks) {
  # get a master set of all nonoverlapping peaks from all peak sets
  all_peaks_gr <- combine_peaks(peaks, method = method, min_overlap = min_overlap)
  
  } else {
    all_peaks_gr <- peaks[[1]]
  }
  all_peaks_df <- as.data.frame(all_peaks_gr) %>%
    dplyr::select(1:5)
  
  # build a logical overlap table indicating which peaks were detected in which samples
  for (i in 1:length(peaks)) {
    sample_name <- names(peaks)[i]
    i_peaks <- peaks[[i]]
    
    all_peaks_df[,sample_name] <- FALSE
    overlaps <- findOverlaps(all_peaks_gr, i_peaks, minoverlap = min_overlap)@from
    all_peaks_df[overlaps,sample_name] <- TRUE
  }
  return(all_peaks_df)
}


# function to filter ChIP peaks to include only those peaks detected in multiple replicates
# peaks should be a GRangesList object
# method should be either overlap or merge
filter_by_replicates <- function(peaks, method = "overlap", min_overlap = 1L) {
  require(tidyverse)
  require(GenomicRanges)
  
  # check input
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GrangesList")
  }

  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  # set names
  if(is.null(names(peaks))) {
    names(peaks) <- paste0("sample_", seq((peaks)))
  }

  # get table of peak overlaps
  overlap_table <- peak_overlap_table(peaks, method = method, min_overlap = min_overlap)
  keep_cols <- names(peaks)
  
  # filter_peaks to include only those detected in all replicates
  filtered_peaks <- overlap_table %>%
    dplyr::mutate(keep = rowSums(select(., all_of(keep_cols)))) %>%
    dplyr::filter(keep == length(peaks)) %>%
    dplyr::select(-all_of(keep_cols), -keep) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
    return(filtered_peaks)
  
}


# function to plot the cumulative distribution function for a set of peaks
# takes as input a Grangeslist object
plot_peak_ECDF <- function(peaks, return_plot = FALSE) {
  require(tidyverse)

  # plot ECDF function
  p <- peaks %>%
    ggplot(aes(log2(fold_enrichment), color = sample)) + 
    stat_ecdf() 
  
  print(p)
  
  if (return_plot) {return(p)}
}

# functions for plotting peak overlap network ----------------------------------
# function to build adjacency matrix of peaks overlap
# peaks should be a GrangesList object containing multiple peak sets to compare
# min_overlap should be an integer specifying the minimum overlap in bp to consider two peaks as overlapping
# type should be 'n' or 'percent' and specify whether the values in the adjacency matrix should be number of peaks or percentages
peak_adj_matrix <- function(peaks, mat_type = "n", ...) {
  # check arguments
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a GRangesList")
  }
  
  
  if (!any(mat_type == c("n", "percent"))) {
    stop("Invalid mat_type provided. mat_type should be 'n' or 'percent' ")
  }
  
  # build adjacency matrix
  overlap_table_n <- matrix(0, ncol=length(peaks), nrow = length(peaks))
  colnames(overlap_table_n) <- names(peaks)
  rownames(overlap_table_n) <- names(peaks)
  overlap_table_p <- overlap_table_n
  
  # compute overlap between all combinations of peak sets and store overlap in tables
  for (i in seq(peaks)) {
    for (j in seq(peaks)) {
      i.gr <- peaks[[i]]
      j.gr <- peaks[[j]]
      
      n_overlapping <- length(subsetByOverlaps(i.gr, j.gr, ...))
      p_overlapping <- length(subsetByOverlaps(i.gr, j.gr, ...)) / length(i.gr)
      
      overlap_table_n[i,j] <- n_overlapping
      overlap_table_p[i,j] <- p_overlapping
    }
    
    
  }
  # return adjacency matrix
  if (mat_type == "n") {
    return(overlap_table_n)
  }
  
  if (mat_type == "percent") {
    return(overlap_table_p)
  }
}
