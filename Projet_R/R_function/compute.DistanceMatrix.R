library(dtwclust)
library(proxy)

compute.DistanceMatrixDTW <-
  function(fenetresViable, normalize = TRUE) {
    matriceDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "dtw_basic",
        upper = FALSE,
        diag = FALSE
      )
    if (normalize) {
      miniDTW <- min(matriceDTW)
      maxiDTW <- max(matriceDTW)
      
      matriceDTW <- matriceDTW - miniDTW
      matriceDTW <- matriceDTW / (maxiDTW - miniDTW)
      
      return(list(matriceDTW, miniDTW, maxiDTW))
    } else{
      return(matriceDTW)
    }
  }

compute.DistanceMatrixSDTW <-
  function(fenetreViable, g = 0.01, normalize = TRUE) {
    matriceSDTW <-
      proxy::dist(
        t(fenetreViable),
        method = "sdtw",
        gamma = g,
        upper = FALSE,
        diag = TRUE
      )
    if (normalize) {
      miniSDTW <- min(matriceSDTW)
      maxiSDTW <- max(matriceSDTW)
      
      matriceSDTW <- matriceSDTW - miniSDTW
      matriceSDTW <- matriceSDTW / (maxiSDTW - miniSDTW)
      
      return(list(matriceSDTW, miniSDTW, maxiSDTW))
    } else{
      return(matriceDTW)
    }
  }