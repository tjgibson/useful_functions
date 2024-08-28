# create custom motif class ====================================================
modisco_motif <- setClass("modisco_motif",
                           
                           slots = list(
                             name = "character",
                             n_seqlets = "numeric",
                             consensus = "character",
                             PPM = "matrix",
                             CWM = "matrix"
                           )
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



# internal function to parse motif =============================================
h5_to_motif <- function(modisco_pattern_list, trim = TRUE, IC_thresh = 0.25) {
  require(stringr)
  require(universalmotif)
  
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
    
    # trim motif based on information content
    if (trim) {
      icm <- convert_type(UM_ppm, "ICM")
      
      IC_pass <- colSums(icm@motif) > IC_thresh
      pass_rle <- rle(IC_pass)
      keep_start <- ifelse(!pass_rle$values[1], pass_rle$lengths[1] + 1 , 1)
      n_runs <- length(pass_rle$lengths)
      keep_end <- ifelse(!pass_rle$values[n_runs], sum(pass_rle$lengths) - pass_rle$lengths[n_runs] , sum(pass_rle$lengths))
      
      ppm <- ppm[,keep_start:keep_end]
      cwm <- cwm[,keep_start:keep_end]
    }
    
    
    # create motif object
    motif <- modisco_motif(
      name = m, 
      n_seqlets = as.numeric(h5_motif$seqlets$n_seqlets), 
      consensus = paste0(colnames(ppm), collapse = ""),
      PPM = ppm, 
      CWM = cwm
    )
    
    motifs[[m]] <- motif
    
  }
  return(motifs)
  }

# parse modisco results file into list of motifs ===============================
read_modisco_h5 <- function(file, neg_patterns = TRUE, trim = TRUE, IC_thresh = 0.1) {
  require(rhdf5)
  require(tidyverse)
  
  # check if input file exists
  if (!file.exists(file)) stop("input file does not exist")
  
  # initiate motif list object
  motifs <- list()
  
  # read motifs with positive contribution scores
  pos_patterns <- h5read(file, "/pos_patterns")
  names(pos_patterns) <- paste0("pos_", names(pos_patterns))
  pos_motifs <- h5_to_motif(pos_patterns, trim = trim, IC_thresh = IC_thresh)
  
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

# function to get reverse complement of motif
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
