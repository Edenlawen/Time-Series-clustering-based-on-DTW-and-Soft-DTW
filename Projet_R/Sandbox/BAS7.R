library(TSA)
library(imputeTS)
library(DTWBI)
library(dtw)
library(ggplot2)
library(readr)
library(gridExtra)
library(Rcpp)
library(cluster)
library(factoextra)
library(tidyr)
library(dplyr)
library(kneedle)
library(dtwclust)
source("R_function/EC_completion.r")
source("R_function/EC_compareCourbe.r")
source("R_function/deformation.r")
source("R_function/soft_DTW_symmetric2.r")

data("google")

#-------------------------------------------------------------------------------
# Etape 1
x <- seq(1:11)
y1 <- c(1, 1, 2, 4, 8, 10, 8, 4, 2, 1, 1)
y2 <- seq(0, 5, by = 0.5)
y3 <- c(1, 1, 1, 1, 2, 4, 8, 10, 8, 4, 2)
y4 <- c(1, 1, 1, 1, 1.5, 10, 1.5, 1, 1, 1, 1)
y5 <- c(1, 1, 2, 4, 8, 8.5, 8, 4, 2, 1, 1) / 5

plot(x, y1, type = "b")
plot(x, y2, type = "b")
plot(x, y3, type = "b")
plot(x, y4, type = "b")
plot(x, y5, type = "b")

#-------------------------------------------------------------------------------
# Etape 2
# DTW

dtw11 <- dtw_basic(y1, y1, norm = "L2")
dtw12 <- dtw_basic(y1, y2, norm = "L2")
dtw13 <- dtw_basic(y1, y3, norm = "L2")
dtw14 <- dtw_basic(y1, y4, norm = "L2")
dtw15 <- dtw_basic(y1, y5, norm = "L2")

dtw22 <- dtw_basic(y2, y2, norm = "L2")
dtw23 <- dtw_basic(y2, y3, norm = "L2")
dtw24 <- dtw_basic(y2, y4, norm = "L2")
dtw25 <- dtw_basic(y2, y5, norm = "L2")

dtw33 <- dtw_basic(y3, y3, norm = "L2")
dtw34 <- dtw_basic(y3, y4, norm = "L2")
dtw35 <- dtw_basic(y3, y5, norm = "L2")

dtw44 <- dtw_basic(y4, y4, norm = "L2")
dtw45 <- dtw_basic(y4, y5, norm = "L2")

dtw55 <- dtw_basic(y5, y5, norm = "L2")

#-------------------------------------------------------------------------------
# Etape 3
# Matrice similiraité
tab <-
  c(dtw11,
    dtw12,
    dtw13,
    dtw14,
    dtw15,
    dtw22,
    dtw23,
    dtw24,
    dtw25,
    dtw33,
    dtw34,
    dtw35,
    dtw44,
    dtw45,
    dtw55)
matrice <- matrix(nrow = 5, ncol = 5)
t <- 1
for (i in 1:5) {
  for (j in i:5) {
    matrice[i, j] <- tab[t]
    matrice[j, i] <- tab[t]
    t <- t + 1
  }
}

print(matrice)
#-------------------------------------------------------------------------------
# Etape 4
# Soft DTW
g <- 0.01

dtw11S <- sdtw(y1, y1, gamma = g)
dtw12S <- sdtw(y1, y2, gamma = g)
dtw13S <- sdtw(y1, y3, gamma = g)
dtw14S <- sdtw(y1, y4, gamma = g)
dtw15S <- sdtw(y1, y5, gamma = g)

dtw22S <- sdtw(y2, y2, gamma = g)
dtw23S <- sdtw(y2, y3, gamma = g)
dtw24S <- sdtw(y2, y4, gamma = g)
dtw25S <- sdtw(y2, y5, gamma = g)

dtw33S <- sdtw(y3, y3, gamma = g)
dtw34S <- sdtw(y3, y4, gamma = g)
dtw35S <- sdtw(y3, y5, gamma = g)

dtw44S <- sdtw(y4, y4, gamma = g)
dtw45S <- sdtw(y4, y5, gamma = g)

dtw55S <- sdtw(y5, y5, gamma = g)

#-------------------------------------------------------------------------------
# Etape 5
# Matrice similarité Soft DTW
tabs <-
  c(dtw11S,
    dtw12S,
    dtw13S,
    dtw14S,
    dtw15S,
    dtw22S,
    dtw23S,
    dtw24S,
    dtw25S,
    dtw33S,
    dtw34S,
    dtw35S,
    dtw44S,
    dtw45S,
    dtw55S)
matriceSoft <- matrix(nrow = 5, ncol = 5)
t <- 1
for (i in 1:5) {
  for (j in i:5) {
    matriceSoft[i, j] <- tabs[t]
    matriceSoft[j, i] <- tabs[t]
    t <- t + 1
  }
}

print(matriceSoft)
#-------------------------------------------------------------------------------
# Autres: Soft DTW test function
g <- 0.0001
result <- sdtw(y1, y2, gamma = g)

#-------------------------------------------------------------------------------
# Autres: Test de la fonction sdtw_cent
data("google")
res_sdtw_cent <- sdtw_cent(google)
df <- data.frame(Base = google, sdtw_res = res_sdtw_cent)

ggplot(df, aes(x = 1:length(google))) +
  geom_line(aes(y = Base, color = "Base")) +
  geom_line(aes(y = sdtw_res, color = "SDTW"))

#-------------------------------------------------------------------------------
# Autres: Test de la fonction cosine pour avoir les query intéressante
data("google")
taille_gap <- 3
start_gap <- 82
cos_conf <- 0.95

AY <-
  gapCreation(google, taille_gap / length(google), start_gap)$output_vector

querys <- seq(4, 100)

df <- data.frame(compute.DTWBI_QF(
  AY,
  acceptedHole = 50,
  smallHole = 0,
  S_Query = querys[1],
  verbose = FALSE
))

df <- setNames(df, paste0("Query = ", querys[1]))

for (i in 2:length(querys)) {
  print(querys[i])
  df <- cbind(
    df,
    compute.DTWBI_QF(
      AY,
      acceptedHole = 50,
      smallHole = 0,
      S_Query = querys[i],
      verbose = FALSE
    )
  )
  colnames(df)[ncol(df)] <- paste0("Query = ", querys[i])
}

dft <- t(df)
dft <- as.data.frame(dft)

res_sdtw_cent <-
  sdtw_cent(dft[, start_gap:(start_gap + taille_gap)], sigma = 0.01)
dfres <- data.frame(Base = google, sdtw_res = google, mean_query_cool = google)
dfres$sdtw_res[start_gap:(start_gap + taille_gap)] <- res_sdtw_cent
dfres$mean <- ts(colMeans(dft))

#-------------------------------------------------------------------------------
# Autres: Test de la fonction Cosine
rescosine <- NULL
queryCool <- NULL
temp1 <- unlist(google[start_gap:(start_gap + taille_gap)])
for (i in 1:length(dft[, 1])) {
  temp2 <- unlist(dft[i, start_gap:(start_gap + taille_gap)])
  rescosine <- c(rescosine, abs(cosine(temp1, temp2)))
  if (abs(cosine(temp1, temp2)) >= cos_conf) {
    queryCool <- c(queryCool, i)
  }
}
print(rescosine)
print(rescosine >= cos_conf)
print(queryCool)

tempsig <- unlist(dft[queryCool[1],start_gap:(start_gap + taille_gap)])
if (length(queryCool) > 1) {
  for (i in 2:length(queryCool)) {
    tempsig <- tempsig + dft[queryCool[i],start_gap:(start_gap + taille_gap)]
  }
  tempsig <- tempsig / 3
  tempsig <- t(tempsig)
}
tempsig <- ts(tempsig)

dfres$mean_query_cool[start_gap:(start_gap + taille_gap)] <- tempsig

#-------------------------------------------------------------------------------
# Affiche Graphe résultat

ggplot(dfres, aes(x = 1:length(google))) +
  geom_line(aes(y = Base, color = "Base"), linewidth = 2) +
  geom_line(aes(y = sdtw_res, color = "SDTW")) +
  geom_line(aes(y = mean, color = "Mean")) +
  geom_line(aes(y = mean_query_cool, color = "Mean query cool")) +
  coord_cartesian(xlim = c(start_gap - 5, start_gap + taille_gap + 5))


#-------------------------------------------------------------------------------
# Affiche matrice
print(matrice)
print(matriceSoft)

#-------------------------------------------------------------------------------
# Affiche Graphe google

ggplot(dfres, aes(x = 1:length(google),y = Base)) +
  geom_line(color = "black") +
  labs(x = "Index", y = "Valeur")

#-------------------------------------------------------------------------------
# Enregistre google dans un csv
dftemp <- data.frame(google)
chemin_fichier <- "google.csv"
write.csv(dftemp, file = chemin_fichier, row.names = FALSE)