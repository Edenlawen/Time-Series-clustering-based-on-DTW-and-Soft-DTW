library(imputeTS)
library(readr)
library(ggplot2)

Pluie <- read_delim(
  "Pluie.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

x <- Pluie
y <- na_interpolation(x, option = "linear")
z <- na_interpolation(x, option = "spline")
v <- na_interpolation(x, option = "stine")

ggplot(Pluie, aes(x = times)) +
  geom_line(aes(y = TSL_mean, color = "Data with gaps"), size = 2) +
  geom_line(aes(y = y$TSL_mean, color = "Linear Interpolation")) +
  geom_point(aes(y = y$TSL_mean, color = "Linear Interpolation")) +
  geom_line(aes(y = z$TSL_mean, color = "Spline Interpolation")) +
  geom_point(aes(y = z$TSL_mean, color = "Spline Interpolation")) +
  geom_line(aes(y = v$TSL_mean, color = "Stine Interpolation")) +
  geom_point(aes(y = v$TSL_mean, color = "Stine Interpolation")) +
  scale_x_continuous(limits = c(Pluie$times[410], Pluie$times[475])) +
  labs(title = "Datas de TSL_mean", x = "Times", y = "Values")

ggplot(Pluie, aes(x = times)) +
  geom_line(aes(y = RSL_mean, color = "Data with gaps"), size = 2) +
  geom_line(aes(y = y$RSL_mean, color = "Linear Interpolation")) +
  geom_point(aes(y = y$RSL_mean, color = "Linear Interpolation")) +
  geom_line(aes(y = z$RSL_mean, color = "Spline Interpolation")) +
  geom_point(aes(y = z$RSL_mean, color = "Spline Interpolation")) +
  geom_line(aes(y = v$RSL_mean, color = "Stine Interpolation")) +
  geom_point(aes(y = v$RSL_mean, color = "Stine Interpolation")) +
  scale_x_continuous(limits = c(Pluie$times[410], Pluie$times[475]))


#Tracé du fichier pluie avec les missing datas
par(mfrow = c(1, 2))

plot(
  Pluie$times,
  Pluie$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'TSL_mean' du fichier Pluie.csv
     \n avec ces missings datas "
)
lines(Pluie$times,
      Pluie$TSL_mean,
      type = "o",
      col = "red")

plot(
  Pluie$times,
  Pluie$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'RSL_mean' du fichier Pluie.csv
     \n avec les missings datas "
)
lines(Pluie$times,
      Pluie$RSL_mean,
      type = "o",
      col = "red")


#Tracé du fichier pluie avec les datas complété avec une interpolation linéaire
par(mfrow = c(1, 2))

plot(
  y$times,
  y$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y$times, y$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

plot(
  y$times,
  y$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y$times, y$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

#Tracé du fichier pluie avec les datas complété avec une interpolation spline
par(mfrow = c(1, 2))

plot(
  z$times,
  z$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(z$times, z$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

plot(
  z$times,
  z$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(z$times, z$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

#Tracé du fichier pluie avec les datas complété avec une interpolation stine
par(mfrow = c(1, 2))

plot(
  v$times,
  v$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(v$times, v$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation stine")

plot(
  v$times,
  v$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(v$times, v$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation stine")


#Tracé du fichier pluie avec les datas complété avec toutes les interpolations
#Mise en place de l'intervale pour regarder les données
a <- 350
b <- 380

par(mfrow = c(2, 2))

plot(
  Pluie$times[a:b],
  Pluie$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'RSL_mean' du fichier Pluie.csv
     \n avec les missings datas "
)
lines(Pluie$times[a:b],
      Pluie$RSL_mean[a:b],
      type = "o",
      col = "red")

plot(
  y$times[a:b],
  y$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y$times[a:b], y$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

plot(
  z$times[a:b],
  z$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(z$times[a:b], z$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

plot(
  v$times[a:b],
  v$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(v$times[a:b], v$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation stine")
