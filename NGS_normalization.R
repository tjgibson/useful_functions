
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
# by default, the function will perform CPM normalization on both bw files prior to computing the ratio. Disable if files are already count normalized.
normalize_bigwigs <- function(input_file, IP_file, psuedocount = 0.01, cpm_normalize = TRUE) {
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
  
  # perform CPM normalization of files
  message("performing CPM normalization")
  if (cpm_normalize) {
    all_intervals$IP_depth <- all_intervals$IP_depth / (sum(all_intervals$IP_depth) / 1e6)
    all_intervals$input_depth <- all_intervals$input_depth / (sum(all_intervals$input_depth) / 1e6)
    
  }
  
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
  all_bins <- tileGenome(seqinfo(bw), tilewidth=min_binsize,cut.last.tile.in.chrom=TRUE)
  
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

# collapse a GRanges object to merge regions with the same score
collapse_gr <- function(gr) {
  collapsed <- unlist(GenomicRanges::reduce(split(gr, ~score)))
  collapsed$score <- as.numeric(names(collapsed))
  names(collapsed) <- NULL
  
  return(collapsed)
  
}
