library(DTWBI)
library(dtw)
library(ggplot2)
source("EC_completion.r")


compareListe <- function(L1, L2, start, taille, verbose = FALSE){
  res <- 0
  for (j in (start):(start+taille)) {
    res <- res + abs(L1[j]-L2[j])
  }
  res <- res / taille
  return(res)
}

searchQuery <- function(Liste, verbose = FALSE){
  query <- 0
  gap <- gapCreation(Liste, rate = 30/4552, begin = 3630)$output_vector
  cond <- 10000
  
  index <- 1:length(Liste)
  
  df <- data.frame(index = index)
  
  df$Ref <- Liste
  
  score <- NULL
  
  for (i in 6:288) {
    df$signal <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                  smallHole = 3, S_Query = i, verbose = FALSE)
    valeur <- compareListe(df$Ref, df$signal, 3630, 30)
    if(verbose==TRUE && i%%20==0){print(i)}
    score <- c(score, valeur)
    if(cond > valeur){
      query = i
      cond = valeur
      if(verbose==TRUE){
        print(i)
        print("TOP")
        print(valeur)
      }
    }
  }
  res <- list(var1 = query, var2 = score)
  return(res)
}

searchQuery2 <- function(Liste, verbose = FALSE){
  query <- 0
  gap <- gapCreation(Liste, rate = 30/4552, begin = 3630)$output_vector
  cond <- 10000
  
  index <- 1:length(Liste)
  
  df <- data.frame(index = index)
  
  df$Ref <- Liste
  
  score <- NULL
  
  for (i in 6:288) {
    df$signal <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                  smallHole = 3, S_Query = i, verbose = FALSE)
    valeur <- compareListe2(df$Ref, df$signal, 3630, 30)
    if(verbose==TRUE && i%%20==0){print(i)}
    score <- c(score, valeur)
    if(cond > valeur){
      query = i
      cond = valeur
      if(verbose==TRUE){
        print(i)
        print("TOP")
        print(valeur)
      }
    }
  }
  res <- list(var1 = query, var2 = score)
  return(res)
}

compareListe2 <- function(L1, L2, start, taille, verbose = FALSE){
  P <- NULL
  Q1 <- NULL
  
  for (j in (start):(start+taille)) {
    P <- c(P,L1[j])
    Q1 <- c(Q1,L2[j])
  }
  
  res <- dtw(x = P, y = Q1, keep.internals = TRUE,  distance.only = FALSE)
  return(res$distance)
}

NH4_Complete <- tsNH4Complete

test <- searchQuery(NH4_Complete,TRUE)
test2 <- searchQuery2(NH4_Complete,TRUE)

df <- data.frame(index = 1:(length(test$var2)))

ggplot(df, aes(x = index)) +
  geom_line(aes(y = test$var2, color = "Complete"))

################################################################################
NH4_Complete <- tsNH4Complete

gap <- gapCreation(tsNH4Complete, rate = 30/4552, begin = 3630)$output_vector

index <- 1:length(NH4_Complete)

df <- data.frame(index = index)

df$NH4_Complete <- NH4_Complete

df$signal2 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                smallHole = 3, S_Query = 2, verbose = FALSE)

df$signal51 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                           smallHole = 3, S_Query = 51, verbose = FALSE)

df$signal60 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                smallHole = 3, S_Query = 60, verbose = FALSE)

ggplot(df, aes(x = index)) +
  geom_line(aes(y = NH4_Complete, color = "Complete"), size = 2) +
  geom_line(aes(y = signal2, color = "Query = 2")) +
  geom_line(aes(y = signal51, color = "Query = 51")) +
  geom_line(aes(y = signal60, color = "Query = 60")) +
  scale_x_continuous(limits = c(index[3625], index[3665])) +
  scale_y_continuous(limits = c(5, 17.5), expand = c(0, 0))

################################################################################
P <- NULL
Q1 <- NULL
Q2 <- NULL

for (i in 3630:3660) {
  P <- c(P,df$NH4_Complete[i])
  Q1 <- c(Q1,df$signal2[i])
  Q2 <- c(Q2,df$signal51[i])
}

comp1 <- dtw(x = P, y = Q1, keep.internals = TRUE,  distance.only = FALSE)
comp2 <- dtw(x = P, y = Q2, keep.internals = TRUE,  distance.only = FALSE)

dtwPlotDensity(comp1)
dtwPlotDensity(comp2)