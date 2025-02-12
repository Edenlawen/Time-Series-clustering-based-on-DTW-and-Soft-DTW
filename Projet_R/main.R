library(TSA)
library(DTWBI)
library(dtw)
library(ggplot2)
library(gridExtra)
library(cluster)
library(factoextra)
library(dtwclust)
library(htmlwidgets)
library(plotly)
library(tidyr)
library(reshape2)
library(tidyverse)
library(reticulate)

rm(list = ls())

source("R_function/EC_completion.r")
source("R_function/EC_compareCourbe.r")
source("R_function/globalF.r")
source("R_function/compute.erreurRMSE.r")
source_python("Python_function/TimeWarp.py")
source("R_function/compute.DistanceMatrix.r")

set.seed(1)


donnee <- datasets::co2
donnee <- unlist(TimeWarp(
  data = r_to_py(donnee),
  nbpts = 500000,
  seq_len = 467,
  condition = 8,
  verbose = FALSE
))

# dataset <- read_csv("csv/bochoiPhuLien.csv")
# dataset <- dataset[,-1]
# dataset <- as.vector(t(dataset))
# dataset <- na.omit(dataset)
# 
# dataset <- r_to_py(dataset)
# 
# rep <-
#   TimeWarp(
#     data = dataset,
#     nbpts = 50000,
#     seq_len = 681,
#     condition = 6,
#     verbose = FALSE
#   )
# donnee <- unlist(rep)

plot(donnee, type = "l")

# Partie Chercher les fenêtres viables pour une taille de query donnée
print("Chercher fenetre viable")
tps <- system.time({
  gapTaille <- 7
  gapStart <- 400
  dataModif <-
    gapCreation(donnee, gapTaille / length(donnee), gapStart)$output_vector
  queryTaille <- 12
  featureRef <-
    globalfeatures(dataModif[(gapStart - queryTaille):(gapStart - 1)])
  queryRef <- dataModif[(gapStart - queryTaille):(gapStart - 1)]
  fenetresViable <- data.frame("queryRef" = queryRef)
  repRef <- donnee[gapStart:(gapStart + gapTaille - 1)]
  reponseViable <- data.frame("repRef" = repRef)
  debut <- 1
  fin <- debut + queryTaille
  
  if (length(donnee) < 1000) {
    borne_a <- 1
    borne_b <- 2
    step_threshold <- sample(borne_a:borne_b, 1)
  } else{
    if (length(donnee) > 10000) {
      borne_a <- 11
      borne_b <- 50
      step_threshold <- sample(borne_a:borne_b, 1)
    } else{
      borne_a <- 3
      borne_b <- 10
      step_threshold <- sample(borne_a:borne_b, 1)
    }
  }
  
  step_threshold <- 2
  threshold_cos <- 0.95
  puiss <- 3
  
  cos_score <- c()
  deb_vect <- c()
  
  while (((fin + queryTaille + gapTaille) < length(donnee))) {
    if (!(
      debut %in% seq(
        gapStart - queryTaille - queryTaille,
        gapStart + gapTaille + queryTaille
      )
    )) {
      featureTemp <- globalfeatures(dataModif[debut:(fin - 1)])
      queryTemp <- dataModif[debut:(fin - 1)]
      cosCompare <- abs(cosine(featureRef, featureTemp))
      
      if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
        cos_score <- c(cos_score, cosCompare)
        deb_vect <- c(deb_vect, debut)
      }
    }
    # step_threshold <- sample(borne_a:borne_b, 1)
    debut <- debut + step_threshold
    fin <- debut + queryTaille
  }
  
  for (i in 1:length(deb_vect)) {
    # print(deb_vect[i])
    debut <- deb_vect[i]
    fin <- debut + queryTaille
    if (compute.indicateurComp(queryRef, dataModif[debut:(fin - 1)])[12] >= 7) {
      fenetresViable <- cbind(fenetresViable,
                              queryTemp <- dataModif[debut:(fin - 1)])
      colnames(fenetresViable)[ncol(fenetresViable)] <-
        paste0("Debut = ", deb_vect[i])
      
      repTemp <- dataModif[fin:(fin + gapTaille - 1)]
      reponseViable <- cbind(reponseViable, repTemp)
      colnames(reponseViable)[ncol(reponseViable)] <-
        paste0("Debut = ", deb_vect[i])
    }
  }
  fenetresViable <- subset(fenetresViable, select = -1)
  reponseViable <- subset(reponseViable, select = -1)
  # rm(debut, fin, cosCompare, featureRef, queryTemp)
})
tps <- tps["elapsed"]
print(paste(
  "Les temps pour la recherche de fenêtre viable est de ",
  tps,
  "secondes"
))

# Partie Calcul de la matrice DTW et Soft-DTW entre toutes les fenêtres viables
# Ca va donner les matrices de distances
print("Calcul des matrices")
tps <- system.time({
  g = 0.001
  print(length(fenetresViable))
  
  matriceDTW <-
    compute.DistanceMatrixDTW(fenetresViable, normalize = FALSE)
  matriceSDTW <-
    compute.DistanceMatrixSDTW(fenetresViable, g, normalize = FALSE)
})
tps <- tps["elapsed"]
print(paste("Les temps pour le calcul des matrice est de ", tps, "secondes"))

#-------------------------------------------------------------------------------
# Matrice de dissimilarity
miniDTW <- min(matriceDTW)
maxiDTW <- max(matriceDTW)
matriceDTW <- matriceDTW - miniDTW
matriceDTW <- matriceDTW / (maxiDTW - miniDTW)

miniSDTW <- min(matriceSDTW)
maxiSDTW <- max(matriceSDTW)
matriceSDTW <- matriceSDTW - miniSDTW
matriceSDTW <- matriceSDTW / (maxiSDTW - miniSDTW)

#-------------------------------------------------------------------------------
# Partie PAM (K-medois)
print("Début de la partie PAM")

nbclusterDTW <-
  fviz_nbclust(t(fenetresViable),
               pam,
               diss = matriceDTW,
               method = "silhouette", )
nbclusterDTW <- nbclusterDTW$data$y
nbclusterDTW <- which.max(nbclusterDTW)
nbclusterSDTW <-
  fviz_nbclust(t(fenetresViable),
               pam,
               diss = matriceSDTW,
               method = "silhouette", )
nbclusterSDTW <- nbclusterSDTW$data$y
nbclusterSDTW <- which.max(nbclusterSDTW)

resultatPamDTW <-
  pam(matriceDTW,
      nbclusterDTW,
      diss = TRUE,
      cluster.only = TRUE)
resultatPamSDTW <-
  pam(matriceSDTW,
      nbclusterSDTW,
      diss = TRUE,
      cluster.only = TRUE)

#   Partie 1: Cluster ayant la moyenne DTW/SDTW la plus faible
avgClusterDTW <- rep(0, times = nbclusterDTW)
avgClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 1:length(fenetresViable)) {
  avgClusterDTW[resultatPamDTW[i]] <-
    avgClusterDTW[resultatPamDTW[i]] + matriceDTW[i]
  avgClusterSDTW[resultatPamSDTW[i]] <-
    avgClusterSDTW[resultatPamSDTW[i]] + matriceSDTW[i]
}

for (i in 1:nbclusterDTW) {
  avgClusterDTW[i] <-
    avgClusterDTW[i] / table(resultatPamDTW)[i]
}

for (i in 1:nbclusterSDTW) {
  avgClusterSDTW[i] <-
    avgClusterSDTW[i] / table(resultatPamSDTW)[i]
}

cat("\nPartie 1\n")
print(paste(
  "Pour DTW, le cluster",
  which.min(avgClusterDTW),
  "a la moyenne de coût DTW la plus faible"
))
print(paste(
  "Pour SDTW, le cluster",
  which.min(avgClusterSDTW),
  "a la moyenne de coût DTW la plus faible"
))

#   Partie 2: Cluster le plus représenté
cat("\nPartie 2\n")
print(paste("Pour DTW, le cluster",
            which.max(table(resultatPamDTW)),
            "a le plus de points"))
print(paste("Pour SDTW, le cluster",
            which.max(table(resultatPamSDTW)),
            "a le plus de points"))

#   Partie 3: Cluster le moyenne d'érreur quadratique la plus basse
avgQuadClusterDTW <- rep(0, times = nbclusterDTW)
avgQuadClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 1:length(fenetresViable)) {
  avgQuadClusterDTW[resultatPamDTW[i]] <-
    avgQuadClusterDTW[resultatPamDTW[i]] + rmse(queryRef, fenetresViable[, i])
  avgQuadClusterSDTW[resultatPamSDTW[i]] <-
    avgQuadClusterSDTW[resultatPamSDTW[i]] + rmse(queryRef, fenetresViable[, i])
}

for (i in 1:nbclusterDTW) {
  avgQuadClusterDTW[i] <-
    avgQuadClusterDTW[i] / table(resultatPamDTW)[i]
}

for (i in 1:nbclusterSDTW) {
  avgQuadClusterSDTW[i] <-
    avgQuadClusterSDTW[i] / table(resultatPamSDTW)[i]
}

cat("\nPartie 3\n")
print(paste(
  "Pour DTW, le cluster",
  which.min(avgQuadClusterDTW),
  "a la moyenne RMSE la plus faible"
))
print(paste(
  "Pour SDTW, le cluster",
  which.min(avgQuadClusterSDTW),
  "a la moyenne RMSE la plus faible"
))

# Partie 4:
avgAmpAvgClusterDTW <- rep(0, times = nbclusterDTW)
avgAmpAvgClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 1:length(reponseViable)) {
  if (queryRef[1] >= fenetresViable[1, i]) {
    avgAmpAvgClusterDTW[resultatPamDTW[i]] <-
      avgAmpAvgClusterDTW[resultatPamDTW[i]] + dtw_basic(queryRef,
                                                         fenetresViable[, i] + mean(abs(queryRef - fenetresViable[, i])))
    avgAmpAvgClusterSDTW[resultatPamSDTW[i]] <-
      avgAmpAvgClusterSDTW[resultatPamSDTW[i]] + sdtw(queryRef,
                                                      fenetresViable[, i] + abs(queryRef - fenetresViable[, i]),
                                                      gamma = g)
  } else{
    avgAmpAvgClusterDTW[resultatPamDTW[i]] <-
      avgAmpAvgClusterDTW[resultatPamDTW[i]] + dtw_basic(queryRef,
                                                         fenetresViable[, i] - abs(queryRef - fenetresViable[, i]))
    avgAmpAvgClusterSDTW[resultatPamSDTW[i]] <-
      avgAmpAvgClusterSDTW[resultatPamSDTW[i]] + sdtw(queryRef,
                                                      fenetresViable[, i] - abs(queryRef - fenetresViable[, i]),
                                                      gamma = g)
  }
}

for (i in 1:nbclusterDTW) {
  avgAmpAvgClusterDTW[i] <-
    avgAmpAvgClusterDTW[i] / table(resultatPamDTW)[i]
}

for (i in 1:nbclusterSDTW) {
  avgAmpAvgClusterSDTW[i] <-
    avgAmpAvgClusterSDTW[i] / table(resultatPamSDTW)[i]
}

cat("\nPartie 4\n")
print(
  paste(
    "Pour DTW, le cluster",
    which.min(avgAmpAvgClusterDTW),
    "a la moyenne de coût DTW la plus faible quand on met le tout les points de la query temp à une distance moyenne entre les 2 querys"
  )
)
print(
  paste(
    "Pour SDTW, le cluster",
    which.min(avgAmpAvgClusterSDTW),
    "a la moyenne de coût DTW la plus faible quand on met le tout les points de la query temp à une distance moyenne entre les 2 querys"
  )
)


# Partie 5:

avgAmpClusterDTW <- rep(0, times = nbclusterDTW)
avgAmpClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 1:length(reponseViable)) {
  if (queryRef[1] >= fenetresViable[1, i]) {
    avgAmpClusterDTW[resultatPamDTW[i]] <-
      avgAmpClusterDTW[resultatPamDTW[i]] + dtw_basic(queryRef,
                                                      fenetresViable[, i] + abs(queryRef[1] - fenetresViable[1, i]))
    avgAmpClusterSDTW[resultatPamSDTW[i]] <-
      avgAmpClusterSDTW[resultatPamSDTW[i]] + sdtw(queryRef,
                                                   fenetresViable[, i] + abs(queryRef[1] - fenetresViable[1, i]),
                                                   gamma = g)
  } else{
    avgAmpClusterDTW[resultatPamDTW[i]] <-
      avgAmpClusterDTW[resultatPamDTW[i]] + dtw_basic(queryRef,
                                                      fenetresViable[, i] - abs(queryRef[1] - fenetresViable[1, i]))
    avgAmpClusterSDTW[resultatPamSDTW[i]] <-
      avgAmpClusterSDTW[resultatPamSDTW[i]] + sdtw(queryRef,
                                                   fenetresViable[, i] - abs(queryRef[1] - fenetresViable[1, i]),
                                                   gamma = g)
  }
}

for (i in 1:nbclusterDTW) {
  avgAmpClusterDTW[i] <-
    avgAmpClusterDTW[i] / table(resultatPamDTW)[i]
}

for (i in 1:nbclusterSDTW) {
  avgAmpClusterSDTW[i] <-
    avgAmpClusterSDTW[i] / table(resultatPamSDTW)[i]
}

cat("\nPartie 5\n")
print(
  paste(
    "Pour DTW, le cluster",
    which.min(avgAmpClusterDTW),
    "a la moyenne de coût DTW la plus faible quand on met le 1er point de la query temp à niveau de la query de réf"
  )
)
print(
  paste(
    "Pour SDTW, le cluster",
    which.min(avgAmpClusterSDTW),
    "a la moyenne de coût DTW la plus faible quand on met le 1er point de la query temp à niveau de la query de réf"
  )
)

# Partie plot de courbe comme Q1; mediane; Q3
dfPAMDTW <- list()
for (i in 1:nbclusterDTW) {
  dfPAMDTW <-
    append(dfPAMDTW, list(data.frame("index" = 1:queryTaille)), after = length(dfPAMDTW))
}
dfPAMSDTW <- list()
for (i in 1:nbclusterSDTW) {
  dfPAMSDTW <-
    append(dfPAMSDTW, list(data.frame("index" = 1:queryTaille)), after = length(dfPAMSDTW))
}
for (i in 1:length(fenetresViable)) {
  df_selected <- dfPAMDTW[[resultatPamDTW[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfPAMDTW[[resultatPamDTW[i]]] <- df_selected
  df_selected <- dfPAMSDTW[[resultatPamSDTW[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfPAMSDTW[[resultatPamSDTW[i]]] <- df_selected
}

for (i in 1:nbclusterDTW) {
  df_selected <- dfPAMDTW[[i]]
  df_selected <- t(df_selected)
  df_selected <- df_selected[2:nrow(df_selected), ]
  dfPAMDTW[[i]] <- df_selected
}
for (i in 1:nbclusterSDTW) {
  df_selected <- dfPAMSDTW[[i]]
  df_selected <- t(df_selected)
  df_selected <- df_selected[2:nrow(df_selected), ]
  dfPAMSDTW[[i]] <- df_selected
}

summaryDTW <- list()
for (i in 1:nbclusterDTW) {
  summaryDTW <-
    append(summaryDTW, list(data.frame("index" = 1:queryTaille)), after = length(summaryDTW))
}

summarySDTW <- list()
for (i in 1:nbclusterSDTW) {
  summarySDTW <-
    append(summarySDTW, list(data.frame("index" = 1:queryTaille)), after = length(summarySDTW))
}


for (i in 1:nbclusterDTW) {
  Q1 <- NULL
  med <- NULL
  Q3 <- NULL
  df_selected <- summaryDTW[[i]]
  df_clust <- dfPAMDTW[[i]]
  if (length(df_clust) != queryTaille) {
    for (j in 1:queryTaille) {
      Q1 <- c(Q1, quantile(df_clust[, j], 0.25))
      med <- c(med, quantile(df_clust[, j], 0.5))
      Q3 <- c(Q3, quantile(df_clust[, j], 0.75))
    }
    df_selected <- cbind(df_selected, queryRef)
    df_selected <- cbind(df_selected, Q1)
    df_selected <- cbind(df_selected, med)
    df_selected <- cbind(df_selected, Q3)
    summaryDTW[[i]] <- df_selected
  } else{
    Q1 <- df_clust[1:queryTaille]
    med <- df_clust[1:queryTaille]
    Q3 <- df_clust[1:queryTaille]
    df_selected <-
      cbind(df_selected, queryRef)
    df_selected <-
      cbind(df_selected, Q1)
    df_selected <-
      cbind(df_selected, med)
    df_selected <-
      cbind(df_selected, Q3)
    summaryDTW[[i]] <- df_selected
  }
}

for (i in 1:nbclusterSDTW) {
  Q1 <- NULL
  med <- NULL
  Q3 <- NULL
  df_selected <- summarySDTW[[i]]
  df_clust <- dfPAMSDTW[[i]]
  if (length(df_clust) != queryTaille) {
    for (j in 1:queryTaille) {
      Q1 <- c(Q1, quantile(df_clust[, j], 0.25))
      med <-
        c(med, quantile(df_clust[, j], 0.5))
      Q3 <-
        c(Q3, quantile(df_clust[, j], 0.75))
    }
    df_selected <-
      cbind(df_selected, queryRef)
    df_selected <-
      cbind(df_selected, Q1)
    df_selected <-
      cbind(df_selected, med)
    df_selected <-
      cbind(df_selected, Q3)
    summarySDTW[[i]] <- df_selected
  } else{
    Q1 <- df_clust[1:queryTaille]
    med <- df_clust[1:queryTaille]
    Q3 <- df_clust[1:queryTaille]
    df_selected <-
      cbind(df_selected, queryRef)
    df_selected <-
      cbind(df_selected, Q1)
    df_selected <-
      cbind(df_selected, med)
    df_selected <-
      cbind(df_selected, Q3)
    summarySDTW[[i]] <- df_selected
  }
  
}

plotPAMDTW <- list()
for (i in 1:nbclusterDTW) {
  t <- summaryDTW[[i]]
  t <- gather(t, key = variable, value = value, -index)
  p <-
    ggplot(data = t, aes(x = index, y = value, color = variable))
  p <-
    p + geom_line() + labs(x = "Jour",
                           y = "Valeur",
                           # title = paste("DTW, PAM, cluster", i))
                           title = paste("DTW, PAM, cluster"))
  
  # Affichage du graphique avec plotly
  plotPAMDTW[[i]] <- ggplotly(p)
}
plotPAMSDTW <- list()
for (i in 1:nbclusterSDTW) {
  t <- summarySDTW[[i]]
  t <- gather(t, key = variable, value = value, -index)
  p <-
    ggplot(data = t, aes(x = index, y = value, color = variable))
  p <-
    p + geom_line() + labs(x = "Jour",
                           y = "Valeur",
                           # title = paste("Soft-DTW, PAM, cluster", i)
                           title = paste("Soft-DTW, PAM, cluster"))
  
  # Affichage du graphique avec plotly
  plotPAMSDTW[[i]] <- ggplotly(p)
}

print(subplot(plotPAMDTW, nrows = 2))
print(subplot(plotPAMSDTW, nrows = 2))

print("Fin de la partie PAM")


# Partie affichage courbe dans le trou
repC1DTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])
repC1SDTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])

repC4DTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])
repC4SDTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])

repC5DTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])
repC5SDTW <-
  data.frame("repRef" = donnee[gapStart:(gapStart + gapTaille - 1)])

for (i in 1:length(reponseViable)) {
  if (resultatPamDTW[i] == which.min(avgClusterDTW)) {
    repC1DTW <- cbind(repC1DTW, reponseViable[, i])
  }
  if (resultatPamSDTW[i] == which.min(avgClusterSDTW)) {
    repC1SDTW <- cbind(repC1SDTW, reponseViable[, i])
  }
  if (resultatPamDTW[i] == which.min(avgAmpAvgClusterDTW)) {
    repC4DTW <- cbind(repC4DTW, reponseViable[, i])
  }
  if (resultatPamSDTW[i] == which.min(avgAmpAvgClusterSDTW)) {
    repC4SDTW <- cbind(repC4SDTW, reponseViable[, i])
  }
  if (resultatPamDTW[i] == which.min(avgAmpClusterDTW)) {
    repC5DTW <- cbind(repC5DTW, reponseViable[, i])
  }
  if (resultatPamSDTW[i] == which.min(avgAmpClusterSDTW)) {
    repC5SDTW <- cbind(repC5SDTW, reponseViable[, i])
  }
}
repC1DTW <- subset(repC1DTW, select = -1)
repC1SDTW <- subset(repC1SDTW, select = -1)
repC4DTW <- subset(repC4DTW, select = -1)
repC4SDTW <- subset(repC4SDTW, select = -1)
repC5DTW <- subset(repC5DTW, select = -1)
repC5SDTW <- subset(repC5SDTW, select = -1)

repC1DTW <- t(repC1DTW)
repC1SDTW <- t(repC1SDTW)
repC4DTW <- t(repC4DTW)
repC4SDTW <- t(repC4SDTW)
repC5DTW <- t(repC5DTW)
repC5SDTW <- t(repC5SDTW)


medRepC1DTW <- vector("numeric", length = 0)
medRepC1SDTW <- vector("numeric", length = 0)

medRepC4DTW <- vector("numeric", length = 0)
medRepC4SDTW <- vector("numeric", length = 0)

medRepC5DTW <- vector("numeric", length = 0)
medRepC5SDTW <- vector("numeric", length = 0)

avgRepC1DTW <- vector("numeric", length = 0)
avgRepC1SDTW <- vector("numeric", length = 0)

avgRepC4DTW <- vector("numeric", length = 0)
avgRepC4SDTW <- vector("numeric", length = 0)

avgRepC5DTW <- vector("numeric", length = 0)
avgRepC5SDTW <- vector("numeric", length = 0)


for (i in 1:gapTaille) {
  # Median
  medRepC1DTW <- c(medRepC1DTW, quantile(repC1DTW[, i], 0.5))
  medRepC1SDTW <-
    c(medRepC1SDTW, quantile(repC1SDTW[, i], 0.5))
  
  medRepC4DTW <- c(medRepC4DTW, quantile(repC4DTW[, i], 0.5))
  medRepC4SDTW <-
    c(medRepC4SDTW, quantile(repC4SDTW[, i], 0.5))
  
  medRepC5DTW <- c(medRepC5DTW, quantile(repC5DTW[, i], 0.5))
  medRepC5SDTW <-
    c(medRepC5SDTW, quantile(repC5SDTW[, i], 0.5))
  
  # Average
  
  avgRepC1DTW <- c(avgRepC1DTW, mean(repC1DTW[, i]))
  avgRepC1SDTW <-
    c(avgRepC1SDTW, mean(repC1SDTW[, i]))
  
  avgRepC4DTW <- c(avgRepC4DTW, mean(repC4DTW[, i]))
  avgRepC4SDTW <-
    c(avgRepC4SDTW, mean(repC4SDTW[, i]))
  
  avgRepC5DTW <- c(avgRepC5DTW, mean(repC5DTW[, i]))
  avgRepC5SDTW <-
    c(avgRepC5SDTW, mean(repC5SDTW[, i]))
}

# medRepC1DTW <- medRepC1DTW + (donnee[gapStart - 1] - medRepC1DTW[1])
# medRepC1SDTW <-
#   medRepC1SDTW + (donnee[gapStart - 1] - medRepC1SDTW[1])
#
# medRepC4DTW <- medRepC4DTW + (donnee[gapStart - 1] - medRepC4DTW[1])
# medRepC4SDTW <-
#   medRepC4SDTW + (donnee[gapStart - 1] - medRepC4SDTW[1])
#
# medRepC5DTW <- medRepC5DTW + (donnee[gapStart - 1] - medRepC5DTW[1])
# medRepC5SDTW <-
#   medRepC5SDTW + (donnee[gapStart - 1] - medRepC5SDTW[1])

df <-
  data.frame(
    "index" = 1:length(medRepC1DTW),
    "main" = donnee[gapStart:(gapStart + gapTaille - 1)],
    # Median
    "medC1DTW" = medRepC1DTW,
    "medC1SDTW" = medRepC1SDTW,
    "medC4DTW" = medRepC4DTW,
    "medC4SDTW" = medRepC4SDTW,
    "medC5DTW" = medRepC5DTW,
    "medC5SDTW" = medRepC5SDTW,
    # Average
    "avgC1DTW" = avgRepC1DTW,
    "avgC1SDTW" = avgRepC1SDTW,
    "avgC4DTW" = avgRepC4DTW,
    "avgC4SDTW" = avgRepC4SDTW,
    "avgC5DTW" = avgRepC5DTW,
    "avgC5SDTW" = avgRepC5SDTW,
    #DTWBI
    "DTWBI" = compute.DTWBI_QF(dataModif, 20, 0, queryTaille, verbose = FALSE)[gapStart:(gapStart +
                                                                                           gapTaille - 1)]
  )
r <- rep(0, times = 7)

cat("\nDTW entre la réponse de référence et la moyenne du cluster DTW selon critére 1\n")
r[1] <- dtw_basic(df$main, df$medC1DTW)
print(r[1])

cat("\nDTW entre la réponse de référence et la moyenne du cluster SDTW selon critére 1\n")
r[2] <- dtw_basic(df$main, df$medC1SDTW)
print(r[2])

cat("\nDTW entre la réponse de référence et la moyenne du cluster DTW selon critére 4\n")
r[3] <- dtw_basic(df$main, df$medC4DTW)
print(r[3])

cat("\nDTW entre la réponse de référence et la moyenne du cluster SDTW selon critére 4\n")
r[4] <- dtw_basic(df$main, df$medC4SDTW)
print(r[4])

cat("\nDTW entre la réponse de référence et la moyenne du cluster DTW selon critére 5\n")
r[5] <- dtw_basic(df$main, df$medC5DTW)
print(r[5])

cat("\nDTW entre la réponse de référence et la moyenne du cluster SDTW selon critére 5\n")
r[6] <- dtw_basic(df$main, df$medC5SDTW)
print(r[6])

cat("\nDTW entre la réponse de référence et DTWBI\n")
r[7] <- dtw_basic(df$main, df$DTWBI)
print(r[7])

t <- gather(df, key = variable, value = value, -index)
p <-
  ggplot(data = t, aes(x = index, y = value, color = variable))
p <-
  p + geom_line() + labs(x = "Jour",
                         y = "Valeur",
                         # title = paste("DTW, PAM, cluster", i))
                         title = paste("Résultat du bouchage"))

# Affichage du graphique avec plotly
print(ggplotly(p))
