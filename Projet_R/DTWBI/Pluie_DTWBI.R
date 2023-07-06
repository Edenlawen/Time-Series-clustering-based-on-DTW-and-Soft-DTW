library(readr)
library(imputeTS)
source("EC_completion.r")

Pluie <- read_delim(
  "Pluie.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

complete_BTWBI <- compute.DTWBI_QF(
  sig = Pluie$RSL_mean,
  # ATTENTION pour "acceptedHole": Si on met une valeur qui est inférieur
  # au plus grand T_gap ce trouvant dans le jeu de donnée, il n'arrivera pas à
  # compléter. Faudrait peut-être automatisé ça en parcourant la liste dans un
  # premier temps pour connaitre la taille max au lieu de la renseigner
  # nous-même.
  acceptedHole = 16,
  smallHole = 4,
  S_Query = 20
)

x1 <- na_interpolation(Pluie)
x2 <- na_interpolation(Pluie, option = "spline")
x3 <- na_interpolation(Pluie, option = "stine")

# Plot de tout le monde pour RSL_mean
ggplot(Pluie, aes(x = times)) +
  geom_line(aes(y = Pluie$RSL_mean, color = "Base"), size = 2) +
  geom_line(aes(y = complete_BTWBI, color = "BTWBI")) +
  geom_point(aes(y = complete_BTWBI, color = "BTWBI")) +
  geom_line(aes(y = x1$RSL_mean, color = "Linear")) +
  geom_point(aes(y = x1$RSL_mean, color = "Linear")) +
  geom_line(aes(y = x2$RSL_mean, color = "Spline")) +
  geom_point(aes(y = x2$RSL_mean, color = "Spline")) +
  geom_line(aes(y = x3$RSL_mean, color = "Stine")) +
  geom_point(aes(y = x3$RSL_mean, color = "Stine")) +
  labs(title = "RSL_mean data from the file 'Pluie.csv'",
       x = "Times", y = "Values")

# Plot de tout le monde pour RSL_mean sur l'intervalle [410;475]
ggplot(Pluie, aes(x = times)) +
  geom_line(aes(y = Pluie$RSL_mean, color = "Base"), size = 2) +
  geom_line(aes(y = complete_BTWBI, color = "BTWBI")) +
  geom_point(aes(y = complete_BTWBI, color = "BTWBI")) +
  geom_line(aes(y = x1$RSL_mean, color = "Linear")) +
  geom_point(aes(y = x1$RSL_mean, color = "Linear")) +
  geom_line(aes(y = x2$RSL_mean, color = "Spline")) +
  geom_point(aes(y = x2$RSL_mean, color = "Spline")) +
  geom_line(aes(y = x3$RSL_mean, color = "Stine")) +
  geom_point(aes(y = x3$RSL_mean, color = "Stine")) +
  scale_x_continuous(limits = c(Pluie$times[410], Pluie$times[475])) +
  scale_y_continuous(limits = c(-54.75,-52.75), expand = c(0, 0)) +
  labs(title = "RSL_mean data from the file 'Pluie.csv'
       in the interval [410;475]", x = "Times", y = "Values")