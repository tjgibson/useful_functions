# General collection of custom R functions used for processing and analyzing genomics data

# deeptools wrappers --------------------------------------------------------------------------------------------
## These functions call deeptools from the command line, but allow for cleaner organization and documentation
## of input files

# deeptools computeMatrix wrapper
dt_compute_matrix <- function(mode, region_files, score_files, matrix_file,run_command = FALSE, ...) {
  require(tidyverse)
  command <- paste("/Users/tylergibson/opt/miniconda3/envs/dt_env/bin/computeMatrix", mode,
                   "-R", str_c(region_files, collapse = " "),
                   "-S", str_c(score_files, collapse = " "),
                   "-o", matrix_file)
  additional_args <- list(...)
  
  command <- paste(command, str_c(additional_args, collapse = " "))
  
  print(command)
  if (run_command == TRUE) {
    system(command)
  }
}

# deeptools plotHeatmap wrapper
dt_plot_heatmap <- function(matrix, out_file,run_command = FALSE, ...) {
  require(tidyverse) 
  command <- paste("/Users/tylergibson/opt/miniconda3/envs/dt_env/bin/plotHeatmap",
                   "-m", matrix,
                   "-out", out_file)
  
  
  additional_args <- list(...)
  
  command <- paste(command, str_c(additional_args, collapse = " "))
  
  print(command)
  
  if (run_command == TRUE) {
    system(command)
  }
}


# NGS normalization functions -------------------------------------------------------------------------------------
# rpkm normalization
rpkm <- function(count_table, widths) {
  require(tidyverse)
  
  row_names <- tibble(id = rownames(count_table))
  
  
  # calculate rpkm
  kb_widths <- widths / 10^3
  
  column_rpkm <- function(x, widths = kb_widths) {
    cr <- (x / widths) / (sum(x) / 10^6 )
    return(cr)
  }
  
  all_rpkm <- count_table %>%
    map(column_rpkm) %>%
    as_tibble() %>%
    bind_cols(row_names) %>%
    column_to_rownames(var = "id")
  
  return(all_rpkm)
  
}

# define function for IP/input normalized bigwig files
# function takes as input two bigwig files corresponding to input and IP from a ChIP experiment
# the function returns a GRanges object with the IP/input ratio
normalize_bigwigs <- function(input_file, IP_file, psuedocount = 0.01) {
  # set dependencies
  require(GenomicRanges)
  require(rtracklayer)
  
  # import input and IP file as bigwig
  message("reading ChIP input file")
  input_gr <- import(input_file)
  message("reading ChIP IP file")
  IP_gr <- import(IP_file)
  
  
  # combine input and IP ranges and match up intervals
  message("matching bigwig bins between input and IP files")
  all_intervals <- disjoin(c(input_gr, IP_gr))
  
  # add the coverage/score for both input and IP
  overlaps <- findOverlaps(all_intervals, input_gr)
  all_intervals$input_depth[overlaps@from] <- input_gr$score[overlaps@to]
  
  overlaps <- findOverlaps(all_intervals, IP_gr)
  all_intervals$IP_depth[overlaps@from] <- IP_gr$score[overlaps@to]
  
  # compute IP/ input ratio
  message("computing IP / input ratio")
  ratio <- all_intervals %>%
    as.data.frame() %>%
    mutate(ratio = case_when(input_depth + IP_depth == 0 ~ 0,
                             input_depth + IP_depth != 0 ~ IP_depth / (input_depth + psuedocount))) %>%
    select(-c("input_depth", "IP_depth")) %>%
    dplyr::rename(score = ratio) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = seqinfo(all_intervals))
  return(ratio)  
}

# zscore normalize bigwigs
zscore_bw <- function(bw) {
  require(tidyverse)
  require(rtracklayer)
  require(GenomicRanges)
  
  if (typeof(bw) == "character") {
    message("reading bigwig file")
    bw <- import(bw)
  }
  
  # for large regions with the same score, expand into equal sized bins
  message("binning genome")
  min_binsize <- min(width(bw))
  all_bins <- tileGenome(seqinfo(bw), tilewidth=10,cut.last.tile.in.chrom=TRUE)
  
  message("getting scores for all bins")
  # add the coverage/score for both input and IP
  all_bins <- subsetByOverlaps(all_bins, bw)
  overlaps <- findOverlaps(all_bins, bw)
  all_bins$score[overlaps@from] <- bw$score[overlaps@to]
  
  # perform z-score normalization
  message("performing z-score normalization")
  all_bins$zscore <- scale(all_bins$score)[,1]
  all_bins$score <- NULL
  all_bins$score <- all_bins$zscore
  all_bins$zscore <- NULL
  # collapse adjacent bins with same score
  collapsed <- unlist(GenomicRanges::reduce(split(all_bins, ~score)))
  collapsed$score <- as.numeric(names(collapsed))
  names(collapsed) <- NULL
  all_bins <- collapsed
  
  #set seqinfo for z-score normalized version
  seqinfo(all_bins) <- seqinfo(bw)
  
  return(all_bins)
}

#functions for processing peaks (e.g. ChIP-seq, ATAC-seq) ----------------------

# function to plot the cumulative distribution function for a set of peaks
# takes as input a character vector of paths to MACS2 output files files (e.g. macs2_peaks.xls)
# if the elements of the character vector are named, these names will be used to label the plot. If not, the filename will be used.
plot_peak_ECDF <- function(in_files, return_plot = FALSE) {
  require(tidyverse)
  # check input file
  if (any(!file.exists(in_files))) {
    stop("narrowPeak file not found")
  }
  
  # set sample names
  if(is.null(names(in_files))) {
    names(in_files) <- str_replace(basename(in_files), "xls", "")
  }
  
  # import peak files 
  all_peaks <- in_files %>%
    map(read_tsv, comment = "#") %>%
    bind_rows(.id = "sample")
  
  # plot ECDF function
  p <- all_peaks %>%
    ggplot(aes(log2(fold_enrichment), color = sample)) + 
    stat_ecdf() 
  
  print(p)
  
  if (return_plot) {return(p)}
}
