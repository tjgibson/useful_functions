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
    track <- file %>%
      rtracklayer::import(which = which) %>%
      `seqlevelsStyle<-`("UCSC") %>%
      DataTrack(...)
    return(track)
  }

  # for bed file, read file and create gViz annotation track
  if(str_detect(file, ".bed")) {
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


