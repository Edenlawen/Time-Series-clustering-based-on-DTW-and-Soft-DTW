library(dtwclust)
library(dtw)
library(proxy)

compute.DistanceMatrixDTW <-
  function(fenetresViable) {
    matriceDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "dtw_basic",
        upper = FALSE,
        diag = FALSE
      )
    return(matriceDTW)
  }

compute.DistanceMatrixSDTW <-
  function(fenetreViable, g = 0.01) {
    matriceSDTW <-
      proxy::dist(
        t(fenetreViable),
        method = "sdtw",
        gamma = g,
        upper = FALSE,
        diag = TRUE
      )
    return(matriceSDTW)
  }