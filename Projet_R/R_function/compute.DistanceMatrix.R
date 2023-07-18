#' Calculates many other datasets from an original to obtain larger ones and different curves
#' Author Merlin Roudier v17/07/2023
#' @param data The original data converted to ndarray for numpy
#' @param nbpts The number of points to be generated. The program will stop when it exceeds this number
#' @param seq_len The size of the number of points desired when a sequence of points will be generated
#' @param verbose
#' @return similarite

compute.DistanceMatrix<-function(fenetreViable, gamma ,verbose=F){
  g = gamma
  if (verbose){print(length(fenetresViable))}
  
  matriceDTW <-
    proxy::dist(
      t(fenetresViable),
      method = "dtw_basic",
      upper = FALSE,
      diag = FALSE
    )
  matriceSDTW <-
    proxy::dist(
      t(fenetresViable),
      method = "sdtw",
      gamma = g,
      upper = FALSE,
      diag = TRUE
    )
  return(matriceDTW,matriceSDTW)
}