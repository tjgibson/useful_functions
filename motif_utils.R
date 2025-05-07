# create custom motif class ====================================================
# validation function
check_modisco_motif <- function(object) {
  errors <- character()
  # check that PPM and CWM have same dimensions
  ppm_dim <- dim(object@PPM)
  cwm_dim <- dim(object@CWM)
  
  if (!all(ppm_dim == cwm_dim)) {
    msg <- paste("PPM and CWM matrices should have same dimensions")
    errors <- c(errors, msg)
  }
  
  # check that length of motif consensus matches PPM and CWM columns
  consensus_length <- nchar(object@consensus)
  
  if (consensus_length != ncol(object@PPM) | consensus_length != ncol(object@CWM)) {
    msg <- paste("number of characters for motif_consensus must match n columns for PPM and CWM matrices")
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

modisco_motif <- setClass("modisco_motif",
                           
                           slots = list(
                             name = "character",
                             n_seqlets = "numeric",
                             consensus = "character",
                             PPM = "matrix",
                             CWM = "matrix"
                           ), 
                          
                          validity = check_modisco_motif
)

setGeneric("motif_name", function(x) standardGeneric("motif_name"))
setMethod("motif_name", "modisco_motif", function(x) x@name)

setGeneric("n_seqlets", function(x) standardGeneric("n_seqlets"))
setMethod("n_seqlets", "modisco_motif", function(x) x@n_seqlets)

setGeneric("motif_consensus", function(x) standardGeneric("motif_consensus"))
setMethod("motif_consensus", "modisco_motif", function(x) x@consensus)

setGeneric("motif_ppm", function(x) standardGeneric("motif_ppm"))
setMethod("motif_ppm", "modisco_motif", function(x) x@PPM)

setGeneric("motif_cwm", function(x) standardGeneric("motif_cwm"))
setMethod("motif_cwm", "modisco_motif", function(x) x@CWM)

# function to convert motif to universalmotif format ===========================
# converts modisco_motif obbject to universalmotif PPM object
modisco_to_universalmotif <- function(motif) {
  um_ppm <- create_motif(
    motif@PPM, 
    type = "PPM", 
    nsites = motif@n_seqlets, 
    name = motif@name
  )
  
  return(um_ppm)
}

# function to trim motif =======================================================
# removes flanking bases with low information content
trim_modisco_motif <- function(motif, trim_by, threshold = 0.25) {
  require(universalmotif)
  
  # check that input object is modisco motif object
  if (class(motif) != "modisco_motif") {stop("provided object is not a modisco_motif object")}
  
  # check that acceptable value is provided for trim_by argument
  if (!(trim_by %in% c("IC", "contribution"))) {stop("trim_by should be set to 'IC' or. 'contribution'")}
  
  # Get bases above threshold based on information content
  if(trim_by == "IC") {
    icm <- universalmotif::convert_type(motif@PPM, "ICM")
    pass_thresh <- colSums(icm@motif) > threshold
  } 
  
  # Get bases above threshold based on contribution score
  else if (trim_by == "contribution") {
    max_contribution <- apply(abs(motif@CWM), 2, FUN = max)
    fraction_max <- max_contribution / (max(max_contribution))
    pass_thresh <- fraction_max > threshold
  }
  
  # trim motif
  pass_rle <- rle(pass_thresh)
  keep_start <- ifelse(!pass_rle$values[1], pass_rle$lengths[1] + 1 , 1)
  n_runs <- length(pass_rle$lengths)
  keep_end <- ifelse(!pass_rle$values[n_runs], sum(pass_rle$lengths) - pass_rle$lengths[n_runs] , sum(pass_rle$lengths))
  
  ppm <- motif@PPM[,keep_start:keep_end]
  cwm <- motif@CWM[,keep_start:keep_end]
  
  
  # create new modisco motif object with trimmed motif
  trimmed_motif <- modisco_motif(
    name = motif@name, 
    n_seqlets = motif@n_seqlets, 
    consensus = paste0(colnames(ppm), collapse = ""),
    PPM = ppm, 
    CWM = cwm
  )
  
  return(trimmed_motif)
}

# internal function to parse motif =============================================
h5_to_motif <- function(modisco_pattern_list, trim = FALSE,trim_by = "IC", trim_thresh = 0.25) {
  require(stringr)
  
  # extract pattern names and sort in order
  motif_names <- names(modisco_pattern_list) |> 
    str_sort(numeric = TRUE)
  
  # initiate motif list
  motifs <- list()
  
  # iterate over motifs and add each to list
  for (m in motif_names) {
    # extract motif PPM and CWM
    h5_motif <- modisco_pattern_list[[m]]
    UM_ppm <- universalmotif::create_motif(h5_motif$sequence, alphabet = "DNA", type = "PPM")
    ppm <- UM_ppm@motif
    cwm <-  h5_motif$contrib_scores
    colnames(cwm) <- colnames(ppm)
    rownames(cwm) <- rownames(ppm)
    
    
    # create motif object
    motif <- modisco_motif(
      name = m, 
      n_seqlets = as.numeric(h5_motif$seqlets$n_seqlets), 
      consensus = paste0(colnames(ppm), collapse = ""),
      PPM = ppm, 
      CWM = cwm
    )
    
    if (trim) {
      motif <- trim_modisco_motif(motif, trim_by = trim_by, threshold = trim_thresh)
    }
    
    motifs[[m]] <- motif
    
  }
  return(motifs)
  }

# parse modisco results file into list of motifs ===============================
read_modisco_h5 <- function(file, neg_patterns = TRUE, trim = FALSE,trim_by = "IC", trim_thresh = 0.25) {
  require(rhdf5)
  require(tidyverse)
  
  # check if input file exists
  if (!file.exists(file)) stop("input file does not exist")
  
  # initiate motif list object
  motifs <- list()
  
  # read motifs with positive contribution scores
  pos_patterns <- h5read(file, "/pos_patterns")
  names(pos_patterns) <- paste0("pos_", names(pos_patterns))
  pos_motifs <- h5_to_motif(pos_patterns, trim = trim, trim_by = trim_by, trim_thresh = trim_thresh)
  
  # read motifs with negative contribution scores
  if (neg_patterns) {
    
    if (any(h5ls(file)$group == "/neg_patterns")) {
    neg_patterns <- h5read(file, "/neg_patterns")
    names(neg_patterns) <- paste0("neg_", names(neg_patterns))
    neg_motifs <- h5_to_motif(neg_patterns, trim = trim, IC_thresh = IC_thresh)
    
    all_motifs <- c(pos_motifs, neg_motifs)
    } else {
      warning(paste("modisco motif file contains no motifs with negative contribution scores:", file))
      all_motifs <- pos_motifs
    }
  } else {
    all_motifs <- pos_motifs
  }
  
  return(all_motifs)
}

# function to convert motif list into dataframe ================================
modisco_list_to_df <- function(modisco_list) {
  require(tidyverse)
  
  # check input data type
  list_class <- class(modisco_list[[1]])
  
  if (list_class != "modisco_motif") {stop("provided object is not a list of modisco_motif objects")}
  
  # initiate df object
  motif_df <- enframe(modisco_list, value = "motif")
  motif_df$n_seqlets <-  map(modisco_list, n_seqlets) |> unlist(use.names = FALSE)
  motif_df$consensus <-  map(modisco_list, motif_consensus) |> unlist(use.names = FALSE)

  return(motif_df)
}

# function to get reverse complement of motif ==================================
modisco_motif_rc <- function(motif) {
  require(Biostrings)
  # check that object class is modisco_motif
  if (class(motif) != "modisco_motif") {stop("provided object is not a modisco_motif object")}
  
  # get consensus RC
  consensus_rc <- motif_consensus(motif) |> 
    DNAString() |> 
    reverseComplement() |> 
    as.character()
  
  # get PPM RC
  ppm_rc <- motif@PPM[4:1,ncol(motif@PPM):1]
  colnames(ppm_rc) <- strsplit(consensus_rc, split = "") |> unlist()
  rownames(ppm_rc) <- c("A", "C", "G", "T")
  
  # get CWM RC
  cwm_rc <- motif@CWM[4:1,ncol(motif@CWM):1]
  colnames(cwm_rc) <- strsplit(consensus_rc, split = "") |> unlist()
  rownames(cwm_rc) <- c("A", "C", "G", "T")
  
  # create RC motif object
  motif_rc <- modisco_motif(
    name = paste0(motif@name, "_RC"), 
    n_seqlets = as.numeric(motif@n_seqlets), 
    consensus = consensus_rc,
    PPM = ppm_rc, 
    CWM = cwm_rc
  )
  
  return(motif_rc)
}


# function to create summary table of modisco motifs ===========================
# summary table contains file paths to temp files of logos plots
# this table can be used to generate an HTML output summary using gt
modisco_summary_table <- function(modisco_table, temp_dir = tempdir()) {
  require(dplyr)
  require(ggseqlogo)

  
  # initialize summary table and generate file paths for logos files
  summary_table <- modisco_table |> 
    mutate(ppm_path = paste0(tempfile(tmpdir = temp_dir),paste0("_",name, "_ppm_logos.png"))) |>
    mutate(ppm_rc_path = paste0(tempfile(tmpdir = temp_dir),paste0("_",name, "_ppm_rc_logos.png"))) |>
    mutate(cwm_path = paste0(tempfile(tmpdir = temp_dir),paste0("_",name, "_cwm_logos.png"))) |> 
    mutate(cwm_rc_path = paste0(tempfile(tmpdir = temp_dir),paste0("_",name, "_cwm_rc_logos.png")))
  
  # write temporary files for logos plots
  for (i in seq(nrow(summary_table))) {
    motif <- summary_table$motif[i][[1]]
    motif_rc <- modisco_motif_rc(motif)
    
    motif_length <- ncol(motif@PPM)
    
    png(summary_table$ppm_path[i], width = (250 * motif_length))
    print(
      ggseqlogo(motif@PPM) +
        theme_void() +
        theme(legend.position = "none")
    )
    dev.off()
    
    png(summary_table$ppm_rc_path[i], width = (250 * motif_length))
    print(
      ggseqlogo(motif_rc@PPM) +
        theme_void() +
        theme(legend.position = "none")
    )
    dev.off()
    

    png(summary_table$cwm_path[i], width = (250 * motif_length))
    print(
      ggseqlogo(motif@CWM, method = "custom") +
        theme_void() +
        theme(legend.position = "none")
    )
    dev.off()
    
    png(summary_table$cwm_rc_path[i], width = (250 * motif_length))
    print(
      ggseqlogo(motif_rc@CWM, method = "custom") +
        theme_void() +
        theme(legend.position = "none")
    )
    dev.off()
    
  }
  
  return(summary_table)
}

# get average contribution score for motifs ====================================
# given a granges and a bigwig file, return the average contribution score across all motifs
add_contr_scores <- function(gr, bw, score_colname = "average_contribution_score") {
  # check that all motifs are of equal width
  motif_width_range <- range(width(gr))
  if (motif_width_range[1] != motif_width_range[2]) {stop("All provided motifs must be the same size")}
  
  # define function to extract matrix from bigwig file
  extract_matrix <- function(track, gr, size, ignore_strand) {
    sum <- suppressWarnings(.Call(
      'BWGFile_summary', path.expand(path(track)),
      as.character(seqnames(gr)), ranges(gr), 
      S4Vectors::recycleIntegerArg(size, "size", length(gr)), 
      "mean", as.numeric(NA_real_), PACKAGE='rtracklayer'
    ))
    M <- do.call( rbind, sum )
    if (!ignore_strand) 
      M[as.character(strand(gr))=='-',] <- M[
        as.character(strand(gr))=='-', ncol(M):1]
    return(M)
  }
  
  # extract matrix of contribution scores
  bwf <- BigWigFile(bw)
  score_matrix <- extract_matrix(bwf, gr, max(width(gr)), ignore_strand = FALSE)
  
  # add contribution scores to granges
  mcols(gr)[score_colname] <- rowMeans(score_matrix)
  
  # return modified granges
  return(gr)
}

# write motifs to modified modisco h5 format ===================================
# this function is intended to prepare motifs for finemo
# Motifs are written in a limited h5 format with only the CWM and PPM
# If motifs are of different lengths, motifs can be optionally padded with 0 values to make them equal length for finemo

pad_motif <- function(modisco_list) {
  # check input data type
  list_class <- class(modisco_list[[1]])
  if (list_class != "modisco_motif") {stop("provided object is not a list of modisco_motif objects")}
  
  # get length of longest motif
  longest_motif_length <- modisco_list |> 
    purrr::map(motif_consensus) |> 
    nchar() |> 
    max()
  
  # iterate over motifs and pad shorter motifs with 0 values
  for (i in seq(modisco_list)) {
    i_motif <- modisco_list[[i]]
    motif_length <- nchar(motif_consensus(i_motif))
    
    # get number of bases to pad on each side
    if (motif_length < longest_motif_length) {
      length_diff <- longest_motif_length - motif_length
      if (length_diff %% 2 == 0) {
        n_left <- length_diff / 2
        n_right <- length_diff / 2
      } else {
        
        n_left <- floor(length_diff / 2)
        n_right <- ceiling(length_diff / 2)
      }
      
      # generate columns for padding
      ppm_left_cols <- matrix(
        data = 0.25, 
        nrow = 4, 
        ncol = n_left, 
        dimnames = list(
          c("A", "C", "G", "T"),
          rep("N", times = n_left)
        )
        )
      
      ppm_right_cols <- matrix(
        data = 0.25, 
        nrow = 4, 
        ncol = n_right, 
        dimnames = list(
          c("A", "C", "G", "T"),
          rep("N", times = n_right)
        )
      )
      
      cwm_left_cols <- matrix(
        data = 0, 
        nrow = 4, 
        ncol = n_left, 
        dimnames = list(
          c("A", "C", "G", "T"),
          rep("N", times = n_left)
        )
      )
      
      cwm_right_cols <- matrix(
        data = 0, 
        nrow = 4, 
        ncol = n_right, 
        dimnames = list(
          c("A", "C", "G", "T"),
          rep("N", times = n_right)
        )
      )
      
      # generate padded motif
      padded_motif <- modisco_motif(
        name = i_motif@name, 
        n_seqlets = i_motif@n_seqlets, 
        consensus = paste0(str_dup("N", n_left), i_motif@consensus, str_dup("N", n_right),  collapse = ""),
        PPM = cbind(ppm_left_cols, i_motif@PPM, ppm_right_cols), 
        CWM = cbind(cwm_left_cols, i_motif@CWM, cwm_right_cols)
      )
      
      modisco_list[[i]] <- padded_motif
    }
  }
  
  return(modisco_list)
}

write_modisco_h5 <- function(modisco_list, out_file) {
  require(rhdf5)
  
  # check input data type
  list_class <- class(modisco_list[[1]])
  
  if (list_class != "modisco_motif") {stop("provided object is not a list of modisco_motif objects")}
  
  # initialise output list
  out_list <- list()
  
  # write motifs to list
  for (i in seq(modisco_list)) {
    i_motif <- modisco_list[[i]]
    i_name <- motif_name(i_motif)
    h5_name <- paste0(i_name, "_", (i - 1))
    
    out_list[[h5_name]] <- list()
    out_list[[h5_name]]$sequence <- motif_ppm(i_motif)
    out_list[[h5_name]]$contrib_scores <- motif_cwm(i_motif)
    
    
  }
  
  if (file.exists(out_file)) {stop("output h5 file already exists")}
  
  h5write(out_list, out_file, name = "pos_patterns")
}
