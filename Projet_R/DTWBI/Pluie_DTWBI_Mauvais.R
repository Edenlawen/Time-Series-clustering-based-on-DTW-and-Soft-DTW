library(readr)
library(imputeTS)
source("EC_completion.r")

Pluie <- read_delim(
  "Pluie.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

Pluie_complete_linear <- na_interpolation(Pluie)
Pluie_complete_spline <- na_interpolation(Pluie, option = "spline")
Pluie_complete_stine <- na_interpolation(Pluie, option = "stine")

query_linear <- gapCreation(Pluie_complete_linear$RSL_mean,
                            rate = 0.05, begin = 293)
query_spline <- gapCreation(Pluie_complete_spline$RSL_mean,
                            rate = 0.05, begin = 293)
query_stine <- gapCreation(Pluie_complete_stine$RSL_mean,
                            rate = 0.05, begin = 293)

t_gap <- query_linear$begin_gap
T_gap <- query_linear$gap_size

results_DTWBI_linear <- DTWBI_univariate(query_linear$output_vector,
                                         t_gap, T_gap)
results_DTWBI_spline <- DTWBI_univariate(query_spline$output_vector,
                                         t_gap, T_gap)
results_DTWBI_stine <- DTWBI_univariate(query_stine$output_vector,
                                        t_gap, T_gap)

# Comparaison entre le jeu de donnée troué et celui complété grace à la fonction
# DTWBI_univariate pour RSL_mean
par(mfrow = c(1, 2))

plot(results_DTWBI$input_vector, type = "o", col = "blue", pch = 16)
lines(results_DTWBI$input_vector, type="o", col="red")

plot(results_DTWBI$output_vector, type = "o", col = "blue", pch = 16)
lines(results_DTWBI$output_vector, type="o", col="red")

# Comparaison entre le jeu de donnée non troué et celui complété grace 
#à la fonction  DTWBI_univariate pour RSL_mean
a <- 285
b <- 325

par(mfrow = c(2, 2))

plot(Pluie_complete$RSL_mean[a:b], type = "o", col = "blue", pch = 16)
lines(Pluie_complete$RSL_mean[a:b], type="o", col="red")
title("Interpolation stine pour RSL_mean")

plot(results_DTWBI$output_vector[a:b], type = "o", col = "blue", pch = 16)
lines(results_DTWBI$output_vector[a:b], type="o", col="red")
title("DTWBI pour RSL_mean")

plot(results_DTWBI$input_vector[a:b], type = "o", col = "blue", pch = 16)
lines(results_DTWBI$input_vector[a:b], type="o", col="red")
title("Missing Data pour RSL_mean")


# Idée de départ:Combler tout les trous qui se présentent au fur et à mesure
#
# Problème: Cependant, ca fonctionne pas car il ne faut qu'un seul
# trou sur la série de donnée pour que la fonction "DTWBI_univariate"
# fonctionne.
#
# Potentiel Solution: Il faudrait prendre une série qu'on compléte en amont
# avec "na_interpolation" et qu'on retienne les endroits où c'était des trous
# avant la complétion pour ensuite créer avec "gapCreation" les gap qu'on
# viendrait alors compléter avec "DTWBI_univariate".
#
# Il y aurait une chance que les données soit faussé car comme il chercherait
# des similitudes avec un graph fait à partir d'interpolations, je sais pas
# vraiment si ca serait juste.

liste_trou_TSL <- list()
liste_taille_gap_TSL <- list()
compteur_TSL <- 0


liste_trou_RSL <- list()
liste_taille_gap_RSL <- list()
compteur_RSL <- 0

for (i in 1:length(Pluie$TSL_mean)) {
  if (is.na(Pluie$TSL_mean[i])) {
    compteur_TSL <- compteur_TSL + 1
    
    if ((is.na(Pluie$TSL_mean[i - 1]) == FALSE && i != 0)
        || (is.na(Pluie$TSL_mean[i]) == TRUE && i == 0)) {
      liste_trou_TSL <- append(liste_trou_TSL, i)
    }
  }
  
  if (is.na(Pluie$RSL_mean[i])) {
    compteur_RSL <- compteur_RSL + 1
    
    if ((is.na(Pluie$RSL_mean[i - 1]) == FALSE && i != 0)
        || (is.na(Pluie$RSL_mean[i]) == TRUE && i == 0)) {
      liste_trou_RSL <- append(liste_trou_RSL, i)
    }
  }
  
  if (i != 0 && i != 1) {
    if (is.na(Pluie$TSL_mean[i - 1]) == TRUE &&
        is.na(Pluie$TSL_mean[i]) == FALSE) {
      temp <- i - liste_trou_TSL[[length(liste_trou_TSL)]]
      liste_taille_gap_TSL <- append(liste_taille_gap_TSL, temp)
    }
    if (is.na(Pluie$RSL_mean[i - 1]) == TRUE &&
        is.na(Pluie$RSL_mean[i]) == FALSE) {
      temp <- i - liste_trou_RSL[[length(liste_trou_RSL)]]
      liste_taille_gap_RSL <- append(liste_taille_gap_RSL, temp)
    }
  }
}

results_DTWBI_TSL <-
  DTWBI_univariate(Pluie$TSL_mean, liste_trou_TSL[[1]],
                   liste_taille_gap_TSL[[1]])
results_DTWBI_RSL <-
  DTWBI_univariate(Pluie$RSL_mean, liste_trou_RSL[[1]],
                   liste_taille_gap_RSL[[1]])