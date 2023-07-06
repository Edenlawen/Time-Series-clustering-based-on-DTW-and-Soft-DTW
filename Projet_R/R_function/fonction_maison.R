#-------------------------------------------------------------------------------
#' @param L1 Vector with signal
#' @param L2 Vector with signal
#' @param start The position where the gap start
#' @param taille Size of the gap
compareListe <- function(L1, L2, start, taille, verbose = FALSE) {
  res <- 0
  for (j in (start):(start + taille)) {
    res <- res + abs(L1[j] - L2[j])
  }
  res <- res / taille
  return(res)
}

#-------------------------------------------------------------------------------
#' @param L1 Vector with signal
#' @param L2 Vector with signal
#' @param start The position where the gap start
#' @param taille Size of the gap

compareListe2 <- function(L1, L2, start, taille, verbose = FALSE) {
  P <- NULL
  Q1 <- NULL
  
  for (j in (start):(start + taille)) {
    P <- c(P, L1[j])
    Q1 <- c(Q1, L2[j])
  }
  
  res <-
    dtw(
      x = P,
      y = Q1,
      keep.internals = TRUE,
      distance.only = FALSE
    )
  return(res$distance)
}

#-------------------------------------------------------------------------------
#' Find the best query
#' @param Liste
#' @param methode Method to find Query: "NAE"; "DTW"
#' @param taille_gap Size of the gap
#' @param start_gap The position where the gap start

searchQuery <- function(Liste,
                        methode = "NAE",
                        taille_gap,
                        start_gap,
                        verbose = FALSE) {
  query <- 0
  gap <-
    gapCreation(Liste, rate = taille_gap / length(Liste), begin = start_gap)$output_vector
  cond <- 10000
  
  index <- 1:length(Liste)
  
  df <- data.frame(index = index)
  
  df$Ref <- Liste
  
  score <- NULL
  
  for (i in 4:60) {
    df$signal <- compute.DTWBI_QF(
      sig = gap,
      acceptedHole = 7,
      smallHole = 0,
      S_Query = i,
      verbose = FALSE
    )
    if (methode == "NAE") {
      valeur <- compareListe(df$Ref, df$signal, start_gap, taille_gap)
    }
    else if (methode == "DTW") {
      valeur <- compareListe2(df$Ref, df$signal, start_gap, taille_gap)
    }
    if (verbose == TRUE && i %% 20 == 0) {
      print(i)
    }
    score <- c(score, valeur)
    if (cond > valeur) {
      query = i
      cond = valeur
      if (verbose == TRUE) {
        print(i)
        print("TOP")
        print(valeur)
      }
    }
  }
  res <- list(var1 = query, var2 = score)
  return(res)
}