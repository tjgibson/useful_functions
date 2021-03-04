
#functions for processing peaks (e.g. ChIP-seq, ATAC-seq) ----------------------
# function to extend peak summits of Granges object
extend_summits <- function(summits_granges, upstream=100, downstream=100) {
  start(summits_granges) <- start(summits_granges) - abs(upstream)
  end(summits_granges) <- end(summits_granges) + abs(downstream)
  return(summits_granges)
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
    stop("one or more peak set is not a Granges object. peaks must be a list of Granges objects")
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
  return(combined_peaks)
  
}

# function to filter ChIP peaks to include only those peaks detected in multiple replicates
# peaks should be a character vector of file paths or a GrangesList object
# if peaks is a list of file paths, files should be bed files or the xls files produced by MACS2 (not a combination of the two)
# method should be either overlap or merge
# enrichment threshold should only be provided if files are xls files produced by MACS2
filter_peaks <- function(peaks, method = "overlap", min_overlap = 1L, output = "filtered_peaks", extend_summits = FALSE, extend_upstream = 100L, extend_downstream = 100L, enrichment_threshold = NULL) {
  require(tidyverse)
  require(GenomicRanges)
  require(rtracklayer)
  
  # check input
  if (length(peaks) < 2) {
    stop("Only one set of peaks provided. For merging peaks, provide multiple peak sets.")
  }
  
  if (is.character(peaks)) {
    if (any(!file.exists(peaks))) {
      stop("one or more files not found")
    }
    
    if (all(str_detect(peaks, ".bed$"))) {
      input_format = "bed"
    }
    else if (all(str_detect(peaks, ".xls$"))) {
      input_format = "macs_xls"
    } else {
      stop("If 'peaks' is a character vector, it should contain a list of bed files or xls files produced by MACS2. All files must be of one type, not a mixture of the two.")
    }
    
  }
  else if (!inherits(peaks, "GRangesList") ) {
    stop("one or more peak set is not a Granges object. peaks must be a list of Granges objects")
  }
  else {
    input_format = "GrangesList"
  }
  
  if (!any(method == c("overlap", "merge"))) {
    stop("Invalid method provided. method should be 'overlap' or 'merge' ")
  }
  
  if (!any(output == c("filtered_peaks", "overlap_table"))) {
    stop("Invalid 'output' provided. 'output' should be 'filtered_peaks' or 'overlap_table' ")
  }
  
  
  if (!is_integer(min_overlap)) {
    stop("min_overlap must be an integer")
  }
  
  if (!is.null(enrichment_threshold) & (input_format != "macs_xls")) {
    stop("enrichment_threshold argument should only be specified if 'peaks is a list of xls files from MACS2'")
  }
  
  
  # if a list of bed files is provided, import them to a GrangesList
  if (input_format == "bed") {
    # set sample names
    if(is.null(names(peaks))) {
      names(peaks) <- str_replace(basename(peaks), ".bed", "")
    }
    
    # import bed files
    peak_set_list <- peaks %>%
      map(rtracklayer::import) %>%
      GRangesList()
    
    if (extend_summits) {
      peak_set_list <- peak_set_list %>%
        extend_summits(upstream = extend_upstream, downstream = extend_downstream)
    }
  }
  
  # if a list of MACS2 XLS files is provided, import them and filter based on fold_enrichment, if necessary
  if (input_format == "macs_xls") {
    # set sample names
    if(is.null(names(peaks))) {
      names(peaks) <- str_replace(basename(peaks), ".xls", "")
    }
    
    # import macs2 xls files
    peak_set_list <- peaks %>%
      purrr::map(read_tsv, comment = "#")
    
    if (!is.null(enrichment_threshold)) {
      message("filtering peaks by fold enrichment")
      peak_set_list <- peak_set_list %>%
        map(dplyr::filter, fold_enrichment > enrichment_threshold)
    }
    if (extend_summits) {
      message("extending summits")
      peak_set_list <- peak_set_list %>%
        purrr::map(select, 1, abs_summit, name) %>%
        purrr::map(makeGRangesFromDataFrame, start.field = "abs_summit", end.field = "abs_summit", keep.extra.columns = TRUE) %>%
        GenomicRanges::GRangesList() %>%
        extend_summits(upstream = extend_upstream, downstream = extend_downstream)
    }
    else {
      peak_set_list <- peak_set_list %>%
        map(select, 1:3) %>%
        map(makeGRangesFromDataFrame, keep.extra.columns = TRUE) %>%
        GRangesList()
    }
  }
  
  if (input_format == "GrangesList") {
    peak_set_list <- peaks
    
    if (extend_summits) {
      peak_set_list <- peak_set_list %>%
        extend_summits(upstream = extend_upstream, downstream = extend_downstream)
    }
  }
  
  # get a master set of all nonoverlapping peaks from all peak sets
  all_peaks_gr <- combine_peaks(peak_set_list, method = method, min_overlap = min_overlap)
  all_peaks_df <- as.data.frame(all_peaks_gr)
  
  # build a logical overlap table indicating which peaks were detected in which samples
  for (i in 1:length(peak_set_list)) {
    sample_name <- names(peak_set_list)[i]
    i_peaks <- peak_set_list[[i]]
    
    all_peaks_df[,sample_name] <- FALSE
    overlaps <- findOverlaps(all_peaks_gr, i_peaks, minoverlap = min_overlap)@from
    all_peaks_df[overlaps,sample_name] <- TRUE
  }
  
  if (output == "overlap_table") {
    return(all_peaks_df)
  }
  
  # filter_peaks to include only those detected in all replicates
  keep_cols <- names(peak_set_list)
  filtered_peaks <- all_peaks_df %>%
    dplyr::mutate(keep = rowSums(select(., all_of(keep_cols)))) %>%
    dplyr::filter(keep == length(peaks)) %>%
    dplyr::select(-all_of(keep_cols)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  if (output == "filtered_peaks") {
    return(filtered_peaks)
  }
  
}


# function to plot the cumulative distribution function for a set of peaks
# takes as input a character vector of paths to MACS2 output files files (e.g. macs2_peaks.xls)
# if the elements of the character vector are named, these names will be used to label the plot. If not, the filename will be used.
read_macs2 <- function(files) {
  # check input file
  if (any(!file.exists(files))) {
    stop("one or more input files not found")
  }
  
  
  # set sample names
  if(is.null(names(files))) {
    names(files) <- str_replace(basename(files), ".xls", "")
  }
  
  # import peak files 
  all_peaks <- files %>%
    map(read_tsv, comment = "#") %>%
    bind_rows(.id = "sample")
  
  return(all_peaks)
}

plot_peak_ECDF <- function(in_files, return_plot = FALSE) {
  require(tidyverse)
  # check input file
  if (any(!file.exists(in_files))) {
    stop("one or more input files not found")
  }
  
  # import peak files
  all_peaks <- read_macs2(in_files)
  
  # plot ECDF function
  p <- all_peaks %>%
    ggplot(aes(log2(fold_enrichment), color = sample)) + 
    stat_ecdf() 
  
  print(p)
  
  if (return_plot) {return(p)}
}
