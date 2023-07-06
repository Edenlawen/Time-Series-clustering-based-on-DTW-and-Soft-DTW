library(TSA)
library(imputeTS)
library(DTWBI)
library(dtw)
library(ggplot2)
library(imputeTS)
library(readr)
library(gridExtra)
library(kableExtra)
library(Rcpp)
source("R_function/EC_completion.r")
source("R_function/fonction_maison.r")
source("R_function/EC_compareCourbe.r")

data("google")

AY <- gapCreation(google, 4 / length(google), 114)$output_vector

test <-
  compute.DTWBI_QF(
    AY,
    acceptedHole = 8,
    smallHole = 0,
    S_Query = 16,
    verbose = TRUE
  )

sig <- data.frame(index = 1:length(google))

sig$complete <- google
sig$DTWBI <- test

ggplot(sig, aes(x = index)) +
  geom_line(aes(y = complete, colour = "Main")) +
  geom_line(aes(y = DTWBI, colour = "DTWBI")) +
  scale_x_continuous(limits = c(sig$index[280], sig$index[310]))
#-------------------------------------------------------------------------------
#' Fonction permettant de rehausser les 2 courbes passer en paramètre pour que
#' aucune valeurs soient en dessous de 0
#' Pour ce faire, je cherche le minimum sur les courbes dans un premier temps.
#' Puis je rajoutes la valeur absolue du minimum à toute les valeurs de 2 courbes
up <- function(courbe1, courbe2) {
  mini <- min(c(courbe1, courbe2))
  if (mini < 0) {
    mini <- abs(mini)
    for (i in 1:length(courbe1)) {
      courbe1[i] <- courbe1[i] + mini
      courbe2[i] <- courbe2[i] + mini
    }
  }
  res <- list(courbe1, courbe2)
  return(res)
}
#-------------------------------------------------------------------------------
ratio <- up(sig$complete, sig$DTWBI)
sig$complete <- ratio[[1]]
sig$DTWBI <- ratio[[2]]

ggplot(sig, aes(x = index)) +
  geom_line(aes(y = complete, colour = "Main")) +
  geom_line(aes(y = DTWBI, colour = "DTWBI")) +
  scale_x_continuous(limits = c(sig$index[280], sig$index[310]))

print(compute.indicateurComp(sig$complete, sig$DTWBI))
GH <- compute.indicateurComp(sig$complete, sig$DTWBI)

#-------------------------------------------------------------------------------
a <- 2
b <- 90
df <- data.frame(
  Query = a:b,
  R2 = NA,
  FA2 = NA,
  FB = NA,
  FS = NA,
  NMSE = NA,
  NMSE_O = NA,
  MG = NA,
  VG = NA,
  simAire = NA,
  distMax = NA,
  all = NA,
  Nombre_de_conditions_vraies = NA,
  Nombre_de_conditions_fausses = NA
)

for (i in a:b) {
  print(i)
  sig$DTWBI <- compute.DTWBI_QF(
    AY,
    acceptedHole = 8,
    smallHole = 0,
    S_Query = i,
    verbose = FALSE
  )
  temp <- compute.indicateurComp(sig$complete, sig$DTWBI, v = F)
  # print(temp)
  for (j in 1:length(temp) + 1) {
    df[i - (a - 1), j] <- temp[[j - 1]]
  }
}