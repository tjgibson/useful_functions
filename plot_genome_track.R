# function for plotting genome browser tracks ==================================
# this function is a convenient wrapper around gviz functions for plotting Drosophila genome browser tracks
file_to_gViz_track <- function(file, which = NULL, ...) {
  require(rtracklayer)
  require(tidyverse)
  require(Gviz)

  # check input
  if (!is.character(file)) stop("files should be a character vector")

  if (!file.exists(file)) stop("one or more input files not found")

  if (!str_detect(file, ".bed|.bw")) stop("Input should be bed file or bigWig file, correct file extensions not detected")


  
  # for bigWig file, read file and create gViz data track
  if(str_detect(file, ".bw")) {
    # check seq level style of input file. If not UCSC, change style of "which" to match target file
    # if (!is.null(which)) {
    #   bwf <- BigWigFile(file)
    #   seqlevelsStyle(which) <- seqlevelsStyle(seqinfo(bwf))
    # }
    # 
    # import file
    track <- file %>%
      rtracklayer::import(which = which) %>%
      `seqlevelsStyle<-`("UCSC") %>%
      DataTrack(...)
    return(track)
  }

  # for bed file, read file and create gViz annotation track
  if(str_detect(file, ".bed")) {
    # check seq level style of input file. If not UCSC, change style of "which" to match target file
    # if (!is.null(which)) {
    #   head_gr <- file %>%
    #     read_tsv(comment = "#", skip = 1, n_max = 5, col_names = c("chromosome", "start", "end")) %>%
    #     makeGRangesFromDataFrame()
    #   
    #   seqlevelsStyle(which) <- seqlevelsStyle(seqinfo(head_gr))
    # }
    
    # import file
    track <- file %>%
      rtracklayer::import(which = which) %>%
      `seqlevelsStyle<-`("UCSC") %>%
      AnnotationTrack(...)
    return(track)
  }

}

plot_genome_tracks <- function(files, chromosome, start, end, track_names=NULL, colors=NULL, groups = NULL, ...) {
  require(rtracklayer)
  require(Gviz)
  require(tidyverse)

  # check input
  if (!is.character(files)) stop("files should be a character vector")

  if (!all(file.exists(files))) stop("one or more input files not found")



  if (!is.character(chromosome)) stop("chromosome should be a character vector")

  if (!is.numeric(start)) stop("start must be numeric")

  if (!is.numeric(end)) stop("end must be numeric")


  if (!is.null(track_names)) {
    if (!is.character(track_names)) stop("track_names should be a character vector")
    if (length(track_names) != length(files)) stop("the number of track_names provided should be equal to the number of input files")
  }

  if (!is.null(colors)) {
    if (!is.character(colors)) stop("colors should be a character vector")

    if (length(colors) != length(files)) stop("the number of colors provided should be equal to the number of input files")
  }

  if (!is.null(groups)) {
    if (!is.character(groups)) stop("groups should be a character vector with the sample groupings")

    if (length(groups) != length(files)) stop("the length of groups should be equal to the number of input files")
  }


    # set track names
  if (is.null(track_names)) {
    track_names <- basename(files)
  }

  # generate track objects for each bigwig file
  region <- GRanges(chromosome, IRanges(start, end))
  
  # check seqinfo of 
    
  message("importing coverage from bigWig files")

  # tracks <- files %>%
  #   map(rtracklayer::import, which = region) %>%
  #   map(`seqlevelsStyle<-`, "UCSC") %>%
  #   imap(~ DataTrack(.x, genome = "dm6", chromosome = chromosome, start=start, end=end, name = track_names[.y],col = colors[.y], col.histogram = colors[.y], fill.histogram = colors[.y])) %>%
  #   set_names(track_names)

  tracks <- files %>%
     imap(~ file_to_gViz_track(.x, which = region, genome = "dm6", chromosome = chromosome, start=start, end=end, name = track_names[.y],col = colors[.y], col.histogram = colors[.y], fill.histogram = colors[.y])) %>%
    set_names(track_names)


  # get axis limits for tracks
  y_limits <- list()
  for (i in 1:length(tracks)) {
    if(inherits(tracks[[i]], "DataTrack")) {
      y_limits[i] <- tracks[[i]] %>% 
        values() %>%
        range() %>%
        extendrange() %>%
        ceiling() %>%
        list()
    } else {
      y_limits[i] <- list(NULL)
      
    }
      
  }
  

  if (!is.null(groups)) {
   for (i in unique(groups)) {
     group_i <- which(i == groups)
     group_min <- min(unlist(y_limits[group_i]))
     group_max <- max(unlist(y_limits[group_i]))

    for (i in group_i) {
      if(inherits(tracks[[i]], "DataTrack")) {y_limits[i] <- list(c(group_min,group_max)) }
      }
   }


  }

  tracks <- tracks %>%
    map2(seq_along(.),(~ `displayPars<-`(x = .x, value = list(ylim = y_limits[[.y]], yTicksAt = y_limits[[.y]] ) ) ) )

  # get UCSC genome annotation to add to plot
  message("getting gene annotations from UCSC")
  knownGenes <- UcscTrack(genome="dm6", chromosome=chromosome, track="NCBI RefSeq", from=start, to=end,
                          trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                          symbol="name2", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes")

  # make list of tracks to include in plot
  all_tracks <- c(tracks, list(knownGenes))


  # plot tracks
  message("plotting tracks")
  plotTracks(all_tracks, from= start, to = end, chromosome = chromosome, type="histogram", transcriptAnnotation="symbol", collapseTranscripts="longest",
             # showAxis = FALSE,
             #background.title = "transparent", col.title="black", col.axis="black", rotation.title = 0,
             ...)


}

# functions for plotting deeptools-style heatmaps ================================================
# function for quickly computing a coverage matrix from a bigwig file. Core code taken from seqPlots package
coverage_matrix <- function(file,regions, upstream = 1000, downstream = 1000, region_width = 1, bin_size =10) {
  require(rtracklayer)
  require(GenomicRanges)
  
  # check arguments
  if (!is.character(file)) stop("file should be a character vector")
  
  if (!file.exists(file)) stop("one or more input files not found")
  
  # get regions to include in heatmap
  gr <- GenomicRanges::resize(regions, region_width, fix='center')
  gr <- suppressWarnings( GenomicRanges::promoters(gr, upstream, downstream ) )
  gr <- trim(gr)
  hm_bins  <- seq(-upstream, downstream, by=bin_size )
  
  # extract_matrix function taken from seqPlots package
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
  
  # read bigwig file to matrix
  bwf <- BigWigFile(file)
  coverage_matrix <- extract_matrix(bwf, gr, length(hm_bins), ignore_strand = TRUE)
  
  return(coverage_matrix)
}

# function for gaussian smoothing of heatmap. Function taken from heatmaps package. Core gaussian smoothing functionality is from EBImage package
gaussianFilter = function(mat, sigma = c(3,3), output.size = dim(mat)) {
  require(EBImage)
  require(heatmaps)
  makeGaussian2D = function(sigma) {
    dim = sigma*3
    x = -dim[1]:dim[1]
    y = -dim[2]:dim[2]
    x_d = vapply(x, dnorm, numeric(1), sd=sigma[1])
    y_d = vapply(y, dnorm, numeric(1), sd=sigma[2])
    brush = outer(x_d, y_d)
    brush
  }
  message("\nApplying Gaussian blur...")
  brush = makeGaussian2D(sigma)
  mat.new = filter2(mat, brush)
  if (!all(output.size == dim(mat))) {
    mat.new = resize(mat.new, output.size[1], output.size[2])
  }
  mat.new
}

plot_heatmap <- function(bw, regions, 
                         upstream = 1000, downstream = 1000, region_width = 1,
                         row_split = NULL,
                         order_by_samples = NULL ,
                         brewer_pal = "Greens",
                         scale_group = NULL,
                         individual_scales = FALSE, 
                         return_heatmap_list = FALSE,
                         min_percentile = 0.01,
                         max_percentile = 0.98,
                         smooth_heatmap = FALSE,
                         row_order = NULL,
                         ...) {
  require(rtracklayer)
  require(GenomicRanges)
  require(tidyverse)
  require(EnrichedHeatmap)
  require(RColorBrewer)
  
  # check arguments
  if (!is.character(bw)) stop("bw should be a character vector")
  
  if (!all(file.exists(bw))) stop("one or more input files not found")
  
  if (any(!str_detect(bw, ".bw$"))) stop("bw file extension missing from one or more files, check input file formats")
  
  
  if (!is.null(scale_group)) {
    if (!is.numeric(scale_group)) stop("scale_group should be a numeric vector")
    if (length(scale_group) != 1 & length(scale_group) != length(bw)) stop("scale_group should be a single number or a vector of numbers of the same length as bw")
  }
  
  # get sample names
  if (is.null(names(bw))) {
    names(bw) <- str_replace(basename(bw), ".bw", "")
  }
  
 # read in coverage matrices
  message("reading coverage matrix from bigWig files")
  mat_list <-bw %>% 
    map(coverage_matrix, regions = regions, upstream = upstream, downstream = downstream, region_width = region_width) 
  
  # check for missing values 
  any_NA <- mat_list %>%
    map(is.na) %>%
    map(any) %>%
    unlist()

  if (any(any_NA)) {
    warning("NAs detected in coverage matrices. Replacing NAs with 0")

    replace_mat_na <- function(mat) {
      mat[is.na(mat)] <- 0
      return(mat)
    }

  mat_list <- mat_list %>%
    map(replace_mat_na)
    }

  
  # get row order
  if (is.null(row_order)) {
  if(is.null(order_by_samples)) {
    order <- mat_list %>% 
      map(as.data.frame) %>% 
      bind_cols %>% 
      rowMeans() %>% 
      order(decreasing = TRUE)
  } else if (length(order_by_samples) == 1) {
    order <- mat_list[[order_by_samples]] %>% 
      rowMeans() %>% 
      order(decreasing = TRUE)
  } else {
    order <- mat_list[order_by_samples] %>% 
      map(as.data.frame) %>%
      bind_cols() %>%
      rowMeans() %>% 
      order(decreasing = TRUE)
  }
  } else {
    order  <- row_order
  }
  
  # # apply gaussian smoothing
  # if (smooth_heatmap) {
  # message("applying gaussian smoothing to heatmap")
  # mat_list <- mat_list %>%
  #   map(gaussianFilter)
  # }
  #   
  message("converting matrix to enrichedHeatmap normalizedMatrix object")
  mat_list <- mat_list %>%
    map(as.normalizedMatrix, k_upstream = upstream/10, k_downstream = downstream/10, k_target = region_width, extend = c(upstream,downstream),
        smooth = smooth_heatmap,
        keep = c(min_percentile, max_percentile))
  
  # get color scale for heatmaps
  if (!individual_scales) {
  col_min <- mat_list %>%
    map(na.exclude) %>%
    map(min) %>%
    unlist() %>%
    min()
  
  col_max <- mat_list %>%
    map(na.exclude) %>%
    map(max) %>%
    unlist() %>%
    max()
  
  col_fun <- circlize::colorRamp2(seq(col_min, col_max, length.out = 9), brewer.pal(9, brewer_pal))
  
  } 
  
  # construct heatmap_list object
  hm_list <- NULL
  
  # set legend
  if (!is.null(row_split)) {
    n_groups <- seq_along(unique(row_split))
    lgd = Legend(at = sort(unique(row_split)), title = "Group", 
                 type = "lines", legend_gp = gpar(col = n_groups))
  } else {
    n_groups <- 1
    lgd <- NULL
  }
  
  for (i in seq_along(mat_list)) {
    # in individual_scales and/or scale_group is set, determine separate scale for individual heatmaps or groups of heatmaps
    if (individual_scales) {
      col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[i]])), max(na.exclude(mat_list[[i]])), length.out = 9), brewer.pal(9, brewer_pal))
    } 
    
    if (!is.null(scale_group)) {
      if (individual_scales) {
        warning("scale_groups and individual_scales are both set. scale_groups will override individual_scales for determining color scales")
      }
      
      if (length(scale_group == length(mat_list))) {
        col_list <- list()
          
          group_i <- which(scale_group == scale_group[i])
          
          col_min <- mat_list[group_i] %>%
            map(na.exclude) %>%
            map(min) %>%
            unlist() %>%
            min()
          
          col_max <- mat_list[group_i] %>%
            map(na.exclude) %>%
            map(max) %>%
            unlist() %>%
            max()
          
          
          
          col_fun <-  circlize::colorRamp2( seq(col_min, col_max, length.out = 9), brewer.pal(9, brewer_pal))
        
        
      } else {
        col_fun <- circlize::colorRamp2(seq(min(na.exclude(mat_list[[scale_group]])), max(na.exclude(mat_list[[scale_group]])), length.out = 9), brewer.pal(9, brewer_pal))
      }
      
    } 

    
    # add heatmap to list
    
    
    message("creating heatmap for sample ", i)
    hm_list <- hm_list + EnrichedHeatmap(mat_list[[i]],  col = col_fun,
                                         row_split = row_split,
                                         row_title_rot = 0,
                                         pos_line = FALSE, 
                                         row_order = order,
                                         column_title = names(bw)[i], 
                                         name = names(bw)[i],
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = n_groups))),
                                         ...)
    
  }
  
  if (!return_heatmap_list) {
  if (is.null(lgd)) {
  draw(hm_list, merge_legends = TRUE)
  } else {
    draw(hm_list,annotation_legend_list = list(lgd), merge_legend = TRUE, heatmap_legend_side = "bottom", 
         annotation_legend_side = "bottom")
  }
  } else {
    return(hm_list)
  }
  
}

