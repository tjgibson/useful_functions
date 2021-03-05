# functions for reading peak files =============================================
# takes a character vector with paths to one or more bed files
# if input is one path, reads bed file into a Granges object
# if input is multiple paths, reads bed files into GrangesList
read_peaks_bed <- function(files, ...) {
  require(rtracklayer)
  require(tidyverse)
  require(GenomicRanges)
  
  # check input
  if (!is.character(files)) stop("files should be a character vector")
  
  if (!all(file.exists(files))) stop("one or more input files not found")
  
  if (!all(str_detect(files, ".bed$"))) warning(".bed file extension missing from one or more files. Check input file format.")
  
  # set names
  if(is.null(names(files))) {
    names(files) <- str_replace(basename(files), ".bed", "")
  }
  
  # read a single bed file
  if (length(files) == 1) {
    gr <- rtracklayer::import(files, ...)
    return(gr)
  }
  
  # read multiple bed files
  if (length(files) > 1) {
    gr <- files %>%
      map(rtracklayer::import, ...) %>%
      GenomicRanges::GRangesList()
    return(gr)
  }
}

# takes a character vector with paths to one or more xls files produced by MACS
# if input is one path, reads bed file into a Granges object
# if input is multiple paths, reads bed files into GrangesList
read_peaks_xls <- function(files, use_summits = FALSE) {
  require(tidyverse)
  require(GenomicRanges)
  
  # check input
  if (!is.character(files)) stop("files should be a character vector")
  
  if (!all(file.exists(files))) stop("one or more input files not found")
  
  if (!all(str_detect(files, ".xls$"))) warning(".xls file extension missing from one or more files. Check input file format.")
  
  # set names
  if(is.null(names(files))) {
    names(files) <- str_replace(basename(files), ".xls", "")
  }
  
  
  # read xls files
  if (!use_summits) {
    gr <- files %>%
      map(read_tsv, comment = "#") %>% 
      map(makeGRangesFromDataFrame, keep.extra.columns = TRUE) %>% 
      GRangesList()
    
  } else {
    # read xls files using peak summits
    gr <- files %>%
      map(read_tsv, comment = "#") %>% 
      map(select, -c("start", "end")) %>%
      map(makeGRangesFromDataFrame, start.field = "abs_summit", end.field = "abs_summit", keep.extra.columns = TRUE) %>%
      GRangesList()
  }
  # reutrn output file
  if (length(files) == 1) return(unlist(gr, use.names = FALSE))
  
  if (length(files) > 1) return(gr)
  
}

# takes a character vector with paths to one or more narrowPeak files
# if input is one path, reads bed file into a Granges object
# if input is multiple paths, reads bed files into GrangesList
read_peaks_narrowPeak <- function(files) {
  require(tidyverse)
  require(GenomicRanges)
  
  # check input
  if (!is.character(files)) stop("files should be a character vector")
  
  if (!all(file.exists(files))) stop("one or more input files not found")
  
  if (!all(str_detect(files, ".narrowPeak$"))) warning(".narrowPeak file extension missing from one or more files. Check input file format.")
  
  # set names
  if(is.null(names(files))) {
    names(files) <- str_replace(basename(files), ".narrowPeak", "")
  }
  
  
  # read narrowPeak files
    gr <- files %>%
      map(read_tsv, 
          col_names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"), 
          col_types = c("ciicncnnni")) %>% 
      map(makeGRangesFromDataFrame, keep.extra.columns = TRUE) %>% 
      GRangesList()
    
  

  # reutrn output file
  if (length(files) == 1) return(unlist(gr))
  
  if (length(files) > 1) return(gr)
  
}


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
peak_overlap_table <- function(peaks, method = "overlap", min_overlap = 1L) {
  
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
  
  
  # get a master set of all nonoverlapping peaks from all peak sets
  all_peaks_gr <- combine_peaks(peaks, method = method, min_overlap = min_overlap)
  all_peaks_df <- as.data.frame(all_peaks_gr) %>%
    select(1:5)
  
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
