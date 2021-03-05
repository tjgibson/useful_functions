# function for plotting genome browser tracks ==================================
# this function is a convenient wrapper around gviz functions for plotting Drosophila genome browser tracks
plot_genome_tracks <- function(files, chromosome, start, end, track_names=NULL, colors=NULL) {
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
  
  # set track names
  if (is.null(track_names)) {
    track_names <- basename(files)
  }
  
  # generate track objects for each bigwig file
  region <- GRanges(chromosome, IRanges(start, end))
  
  message("importing coverage from bigWig files")
  tracks <- files %>%
    map(rtracklayer::import, which = region) %>%
    map(`seqlevelsStyle<-`, "UCSC") %>%
    imap(~ DataTrack(.x, genome = "dm6", chromosome = chromosome, start=start, end=end, name = track_names[.y], col = colors[.y], fill = colors[.y])) %>%
    set_names(track_names)
  
  
  # get UCSC genome annotation to add to plot
  message("getting gene annotations from UCSC")
  knownGenes <- UcscTrack(genome="dm6", chromosome=chrom, track="NCBI RefSeq", from=start, to=end,
                          trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                          symbol="name2", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes")
  
  # make list of tracks to include in plot
  all_tracks <- c(tracks, list(knownGenes))
  
  
  # plot tracks
  message("plotting tracks")
  plotTracks(all_tracks, from= start, to = end, chromosome = chromosome, type="h", transcriptAnnotation="symbol", collapseTranscripts="longest", background.title = "transparent", col.title="black", col.axis="transparent", rotation.title = 0)
  
  
}


