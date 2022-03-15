## set working directory
setwd("/Volumes/TG1/genomics/2021_Gaskill_Gibson/data/2018_12_MG_CHIPseq/R_analysis/")

## load required packages
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(pheatmap)

# import peak files
in_files <- c(GAF_S2 = "./output_data/hc_S2_GAF_IP_peaks.bed",
              GAF_St3_WT = "./output_data/hc_St3_GAF_IP_GFP_GAF_peaks.bed",
              GAF_St5_WT = "./output_data/hc_St5_GAF_IP_GFP_GAF_peaks.bed",
              GAF_St5_ZldRNAi = "./output_data/hc_St5_GAF_IP_Zld_RNAi_peaks.bed",
              Zld_nc8_WT = "./output_data/hc_zld_nc8_peaks.bed",
              Zld_nc13_WT = "./output_data/hc_zld_nc13_peaks.bed",
              Zld_nc14_WT = "./output_data/hc_zld_nc14_peaks.bed",
              CLAMP_nc14_male = "../../../../published_data/Rieder_2019/R_analysis/output_data/nc14_male_aCLAMP_peaks.bed",
              CLAMP_nc14_female = "../../../../published_data/Rieder_2019/R_analysis/output_data/nc14_female_aCLAMP_peaks.bed")


peaks <- in_files %>%
  map(rtracklayer::import) %>% 
  GRangesList()

adj_m <- peak_adj_matrix(peaks)
diag(adj_m) <- 0

library(igraph)

g <- graph_from_adjacency_matrix(adj_m,mode = "undirected", weighted = TRUE)
V(g)$size = elementNROWS(peaks) / 500

plot(g, 
     edge.arrow.size=.4, 
     edge.width = E(g)$weight / 500)


adj_m <- peak_adj_matrix(peaks, type = "percent")
diag(adj_m) <- 0


g <- graph_from_adjacency_matrix(adj_m,mode = "undirected", weighted = TRUE)
#V(g)$size = elementNROWS(peaks) / 500

plot(g, 
     edge.arrow.size=.4, 
     edge.width = E(g)$weight * 10)
