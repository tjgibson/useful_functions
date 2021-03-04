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

