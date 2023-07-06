library(imputeTS)
library(readr)
library(ggplot2)

NH4 <- read_delim(
  "NH4.csv",
  delim = ";",
  escape_double = FALSE,
  trim_ws = TRUE
)

y1 <- na_interpolation(NH4, option = "linear")
y2 <- na_interpolation(NH4, option = "spline")
y3 <- na_interpolation(NH4, option = "stine")

ggplot(NH4, aes(x = times)) +
  geom_line(aes(y = TSL_mean, color = "Data with gaps"), size = 2) +
  geom_line(aes(y = y1$TSL_mean, color = "Linear Interpolation")) +
  geom_point(aes(y = y1$TSL_mean, color = "Linear Interpolation")) +
  geom_line(aes(y = y2$TSL_mean, color = "Spline Interpolation")) +
  geom_point(aes(y = y2$TSL_mean, color = "Spline Interpolation")) +
  geom_line(aes(y = y3$TSL_mean, color = "Stine Interpolation")) +
  geom_point(aes(y = y3$TSL_mean, color = "Stine Interpolation")) +
  scale_x_continuous(limits = c(NH4$times[150], NH4$times[192])) +
  labs(title = "Datas TSL_mean of the file 'NH4.csv'", x = "Times", y = "Values")

ggplot(NH4, aes(x = times)) +
  geom_line(aes(y = RSL_mean, color = "Data with gaps"), size = 2) +
  geom_line(aes(y = y1$RSL_mean, color = "Linear Interpolation")) +
  geom_point(aes(y = y1$RSL_mean, color = "Linear Interpolation")) +
  geom_line(aes(y = y2$RSL_mean, color = "Spline Interpolation")) +
  geom_point(aes(y = y2$RSL_mean, color = "Spline Interpolation")) +
  geom_line(aes(y = y3$RSL_mean, color = "Stine Interpolation")) +
  geom_point(aes(y = y3$RSL_mean, color = "Stine Interpolation")) +
  scale_x_continuous(limits = c(NH4$times[150], NH4$times[192])) +
  scale_y_continuous(limits = c(-39.25, -37.5), expand = c(0, 0)) +
  labs(title = "Datas RSL_mean of the file 'NH4.csv'", x = "Times", y = "Values")

#Tracé du fichier NH4 avec les missing datas
par(mfrow = c(1, 2))

plot(
  x$times,
  x$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'TSL_mean' du fichier NH4.csv
     \n avec ces missings datas "
)
lines(x$times, x$TSL_mean, type = "o", col = "red")

plot(
  x$times,
  x$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'RSL_mean' du fichier NH4.csv
     \n avec les missings datas "
)
lines(x$times, x$RSL_mean, type = "o", col = "red")

#Tracé du fichier NH4 avec les datas complété avec une interpolation linéaire
par(mfrow = c(1, 2))

plot(
  y1$times,
  y1$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y1$times, y1$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

plot(
  y1$times,
  y1$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y1$times, y1$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

#Tracé du fichier NH4 avec les datas complété avec une interpolation spline
par(mfrow = c(1, 2))

plot(
  y2$times,
  y2$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y2$times, y2$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

plot(
  y2$times,
  y2$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y2$times, y2$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

#Tracé du fichier pluie avec les datas complété avec une interpolation stine
par(mfrow = c(1, 2))

plot(
  y3$times,
  y3$TSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y3$times, y3$TSL_mean, type = "o", col = "red")
title("Plot de 'TSL_mean' avec ces datas compléter
     \n avec une interpolation stine")

plot(
  y3$times,
  y3$RSL_mean,
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y3$times, y3$RSL_mean, type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation stine")

#Tracé du fichier NH4 avec les datas complété avec toutes les interpolations
#Mise en place de l'intervale pour regarder les données
a <- 137
b <- 192

par(mfrow = c(2, 2))

plot(
  x$times[a:b],
  x$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value",
  main = "Plot de 'RSL_mean' du fichier NH4.csv
     \n avec les missings datas "
)
lines(x$times[a:b], x$RSL_mean[a:b], type = "o", col = "red")

plot(
  y1$times[a:b],
  y1$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y1$times[a:b], y1$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation linéaire")

plot(
  y2$times[a:b],
  y2$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y2$times[a:b], y2$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation spline")

plot(
  y3$times[a:b],
  y3$RSL_mean[a:b],
  type = "o",
  col = "blue",
  pch = 16,
  xlab = "Temps",
  ylab = "Value"
)
lines(y3$times[a:b], y3$RSL_mean[a:b], type = "o", col = "red")
title("Plot de 'RSL_mean' avec ces datas compléter
     \n avec une interpolation stine")