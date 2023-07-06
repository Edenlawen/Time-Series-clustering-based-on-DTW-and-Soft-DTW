library(TSA)
library(DTWBI)
library(dtw)
library(ggplot2)
library(gridExtra)
library(cluster)
library(factoextra)
library(kneedle)
library(dtwclust)
library(htmlwidgets)
library(plotly)
library(tidyr)
library(reshape2)
library(tidyverse)

rm(list = ls())

source("R_function/compute.erreurRMSE.r")

tp <- system.time({
  pts <- 10
  x <- seq(0, pi, length.out = 12)
  
  # Calculer les valeurs de la sinusoïde
  ref <- c(rep(0, pts), sin(x) * 5, rep(0, pts))
  
  plot(ref, type = "b")
  
  
  amplitude = 5
  nombre_points = 6
  
  # Générer la séquence linéaire
  sequence <- seq(0, amplitude, length.out = nombre_points)
  
  # Répéter la séquence
  signal_triangle <- c(rep(0, pts), sequence, rev(sequence), rep(0, pts))
  
  plot(signal_triangle, type = "b")
  
  #-------------------------------------------------------------------------------
  df <- data.frame("1" = signal_triangle, "2" = ref)
  
  compt <- 1
  
  for (i in 6:20) {
    for (k in 0:pts) {
      sequence <- seq(0, i*1.5, length.out = nombre_points)
      signal_triangle <-
        c(rep(0, k), sequence, rev(sequence), rep(0, (pts*2) - k))
      df <- cbind(df, signal_triangle)
      colnames(df)[ncol(df)] <-
        paste0("Triangle = ", compt)
      compt <- compt + 1
    }
  }
  
  compt <- 1
  
  for (i in 6:40) {
    for (k in 0:pts) {
      x <- seq(0, pi, length.out = 12)
      y <- c(rep(0, k), sin(x) * i * 1.5, rep(0, (pts*2) - k))
      df <- cbind(df, y)
      colnames(df)[ncol(df)] <-
        paste0("Cercle = ", compt)
      compt <- compt + 1
    }
  }
  
  
  
  fenetresViable <- df
  
  fenetresViable <- subset(fenetresViable, select = -1)
  #-------------------------------------------------------------------------------
  # Partie PAM (K-medois)
  print("Début de la partie PAM")
  
  nbclusterDTW <-
    fviz_nbclust(t(fenetresViable),
                 pam,
                 method = "silhouette",
                 k.max = 50)
  nbclusterDTW <- nbclusterDTW$data$y
  nbclusterDTW <- which.max(nbclusterDTW)
  
  
  resultatPamDTW <- pam(t(fenetresViable), nbclusterDTW, diss = FALSE)
  
  # Réduction de dimension avec PCA
  # pcaDTW <- prcomp(matriceDTWdist)
  # pcaSDTW <- prcomp(matriceSDTWdist)
  #
  # # Utilisation des coordonnées PCA pour visualiser les clusters
  # plot3 <-
  #   fviz_pca_ind(pcaDTW, habillage = as.factor(resultatPamDTW$clustering)) + ggtitle("K-Médoïdes avec DTW")
  # plot4 <-
  #   fviz_pca_ind(pcaSDTW, habillage = as.factor(resultatPamSDTW$clustering)) + ggtitle("K-Médoïdes avec Soft-DTW")
  #
  # rm(pcaDTW, pcaSDTW)
  
  #   Partie 3: Cluster le moyenne d'érreur quadratique la plus basse
  avgQuadClusterDTW <- rep(0, times = nbclusterDTW)
  
  for (i in 1:length(fenetresViable)) {
    avgQuadClusterDTW[resultatPamDTW$clustering[i]] <-
      avgQuadClusterDTW[resultatPamDTW$clustering[i]] + rmse(ref, fenetresViable[, i])
  }
  
  for (i in 1:nbclusterDTW) {
    avgQuadClusterDTW[i] <-
      avgQuadClusterDTW[i] / table(resultatPamDTW$clustering)[i]
  }
  
  print(paste(
    "Pour DTW, le cluster",
    which.min(avgQuadClusterDTW),
    "a la moyenne RMSE la plus faible"
  ))
  
  queryTaille <- length(ref)
  
  # Partie plot de courbe comme Q1; mediane; Q3
  dfPAMDTW <- list()
  for (i in 1:nbclusterDTW) {
    dfPAMDTW <-
      append(dfPAMDTW, list(data.frame("index" = 1:queryTaille)), after = length(dfPAMDTW))
  }
  
  for (i in 1:length(fenetresViable)) {
    df_selected <- dfPAMDTW[[resultatPamDTW$clustering[i]]]
    df_selected <- cbind(df_selected, fenetresViable[i])
    dfPAMDTW[[resultatPamDTW$clustering[i]]] <- df_selected
  }
  
  for (i in 1:nbclusterDTW) {
    df_selected <- dfPAMDTW[[i]]
    df_selected <- t(df_selected)
    dfPAMDTW[[i]] <- df_selected
  }
  
  summaryDTW <- list()
  for (i in 1:nbclusterDTW) {
    summaryDTW <-
      append(summaryDTW, list(data.frame("index" = 1:queryTaille)), after = length(summaryDTW))
  }
  
  
  for (i in 1:nbclusterDTW) {
    Q1 <- NULL
    med <- NULL
    Q3 <- NULL
    # moy <- NULL
    df_selected <- summaryDTW[[i]]
    df_clust <- dfPAMDTW[[i]]
    for (j in 1:(queryTaille)) {
      Q1 <- c(Q1, quantile(df_clust[2:length(df_clust[, j]), j], 0.25))
      med <- c(med, quantile(df_clust[, j], 0.5))
      Q3 <- c(Q3, quantile(df_clust[, j], 0.75))
    }
    df_selected <- cbind(df_selected, ref)
    df_selected <- cbind(df_selected, Q1)
    df_selected <- cbind(df_selected, med)
    df_selected <- cbind(df_selected, Q3)
    summaryDTW[[i]] <- df_selected
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
                             title = paste("DTW, PAM, cluster", i))
    plotPAMDTW[[i]] <- ggplotly(p)
  }
  
  plot1 <- subplot(plotPAMDTW, nrows = 3)
  
  print("Fin de la partie PAM SANS DTW")
  
  #-------------------------------------------------------------------------------
  # Partie Calcul de la matrice DTW et Soft-DTW entre toutes les fenêtres viables
  # Ca va donner les matrices de distances
  print("Calcul des matrices")
  tps <- system.time({
    g = 0.001
    print(length(fenetresViable))
    
    matriceDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "dtw_basic",
        upper = FALSE,
        diag = FALSE
      )
    matriceSDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "sdtw",
        gamma = g,
        upper = FALSE,
        diag = FALSE
      )
  })
  tps <- tps["elapsed"]
  print(paste("Les temps pour le calcul des matrice est de ", tps, "secondes"))
  
  
  #-------------------------------------------------------------------------------
  # Normalisation
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
    fviz_nbclust(
      t(fenetresViable),
      pam,
      diss = matriceDTW,
      method = "silhouette",
      k.max = 25
    )
  nbclusterDTW <- nbclusterDTW$data$y
  nbclusterDTW <- which.max(nbclusterDTW)
  nbclusterSDTW <-
    fviz_nbclust(
      t(fenetresViable),
      pam,
      diss = matriceSDTW,
      method = "silhouette",
      k.max = 25
    )
  nbclusterSDTW <- nbclusterSDTW$data$y
  nbclusterSDTW <- which.max(nbclusterSDTW)
  
  # nbclusterSDTW <- nbclusterDTW
  
  resultatPamDTW <- pam(matriceDTW, nbclusterDTW, diss = TRUE)
  resultatPamSDTW <- pam(matriceSDTW, nbclusterSDTW, diss = TRUE)
  
  # Réduction de dimension avec PCA
  # pcaDTW <- prcomp(matriceDTWdist)
  # pcaSDTW <- prcomp(matriceSDTWdist)
  #
  # # Utilisation des coordonnées PCA pour visualiser les clusters
  # plot3 <-
  #   fviz_pca_ind(pcaDTW, habillage = as.factor(resultatPamDTW$clustering)) + ggtitle("K-Médoïdes avec DTW")
  # plot4 <-
  #   fviz_pca_ind(pcaSDTW, habillage = as.factor(resultatPamSDTW$clustering)) + ggtitle("K-Médoïdes avec Soft-DTW")
  #
  # rm(pcaDTW, pcaSDTW)
  
  #   Partie 1: Cluster ayant la moyenne DTW/SDTW la plus faible
  avgClusterDTW <- rep(0, times = nbclusterDTW)
  avgClusterSDTW <- rep(0, times = nbclusterSDTW)
  
  for (i in 1:length(fenetresViable)) {
    avgClusterDTW[resultatPamDTW$clustering[i]] <-
      avgClusterDTW[resultatPamDTW$clustering[i]] + matriceDTW[i]
    avgClusterSDTW[resultatPamSDTW$clustering[i]] <-
      avgClusterSDTW[resultatPamSDTW$clustering[i]] + matriceSDTW[i]
  }
  
  for (i in 1:nbclusterDTW) {
    avgClusterDTW[i] <-
      avgClusterDTW[i] / table(resultatPamDTW$clustering)[i]
  }
  
  for (i in 1:nbclusterSDTW) {
    avgClusterSDTW[i] <-
      avgClusterSDTW[i] / table(resultatPamSDTW$clustering)[i]
  }
  
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
  print(paste("Pour DTW, le cluster", which.max(table(
    resultatPamDTW$clustering
  )), "a le plus de points"))
  print(paste("Pour SDTW, le cluster", which.max(table(
    resultatPamSDTW$clustering
  )), "a le plus de points"))
  
  #   Partie 3: Cluster le moyenne d'érreur quadratique la plus basse
  avgQuadClusterDTW <- rep(0, times = nbclusterDTW)
  avgQuadClusterSDTW <- rep(0, times = nbclusterSDTW)
  
  for (i in 1:length(fenetresViable)) {
    avgQuadClusterDTW[resultatPamDTW$clustering[i]] <-
      avgQuadClusterDTW[resultatPamDTW$clustering[i]] + rmse(ref, fenetresViable[, i])
    avgQuadClusterSDTW[resultatPamSDTW$clustering[i]] <-
      avgQuadClusterSDTW[resultatPamSDTW$clustering[i]] + rmse(ref, fenetresViable[, i])
  }
  
  for (i in 1:nbclusterDTW) {
    avgQuadClusterDTW[i] <-
      avgQuadClusterDTW[i] / table(resultatPamDTW$clustering)[i]
  }
  
  for (i in 1:nbclusterSDTW) {
    avgQuadClusterSDTW[i] <-
      avgQuadClusterSDTW[i] / table(resultatPamSDTW$clustering)[i]
  }
  
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
    df_selected <- dfPAMDTW[[resultatPamDTW$clustering[i]]]
    df_selected <- cbind(df_selected, fenetresViable[i])
    dfPAMDTW[[resultatPamDTW$clustering[i]]] <- df_selected
    df_selected <- dfPAMSDTW[[resultatPamSDTW$clustering[i]]]
    df_selected <- cbind(df_selected, fenetresViable[i])
    dfPAMSDTW[[resultatPamSDTW$clustering[i]]] <- df_selected
  }
  
  for (i in 1:nbclusterDTW) {
    df_selected <- dfPAMDTW[[i]]
    df_selected <- t(df_selected)
    df_selected <- df_selected[2:nrow(df_selected),]
    dfPAMDTW[[i]] <- df_selected
  }
  for (i in 1:nbclusterSDTW) {
    df_selected <- dfPAMSDTW[[i]]
    df_selected <- t(df_selected)
    df_selected <- df_selected[2:nrow(df_selected),]
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
    moy <- NULL
    df_selected <- summaryDTW[[i]]
    df_clust <- dfPAMDTW[[i]]
    for (j in 1:queryTaille) {
      Q1 <- c(Q1, quantile(df_clust[, j], 0.25))
      med <- c(med, quantile(df_clust[, j], 0.5))
      Q3 <- c(Q3, quantile(df_clust[, j], 0.75))
      moy <- c(moy, mean(df_clust[, j]))
    }
    df_selected <- cbind(df_selected, ref)
    df_selected <- cbind(df_selected, Q1)
    df_selected <- cbind(df_selected, med)
    df_selected <- cbind(df_selected, Q3)
    df_selected <- cbind(df_selected, moy)
    summaryDTW[[i]] <- df_selected
  }
  
  for (i in 1:nbclusterSDTW) {
    Q1 <- NULL
    med <- NULL
    Q3 <- NULL
    moy <- NULL
    df_selected <- summarySDTW[[i]]
    df_clust <- dfPAMSDTW[[i]]
    for (j in 1:queryTaille) {
      Q1 <- c(Q1, quantile(df_clust[, j], 0.25))
      med <- c(med, quantile(df_clust[, j], 0.5))
      Q3 <- c(Q3, quantile(df_clust[, j], 0.75))
      moy <- c(moy, mean(df_clust[, j]))
    }
    df_selected <- cbind(df_selected, ref)
    df_selected <- cbind(df_selected, Q1)
    df_selected <- cbind(df_selected, med)
    df_selected <- cbind(df_selected, Q3)
    df_selected <- cbind(df_selected, moy)
    summarySDTW[[i]] <- df_selected
  }
  
  plotPAMDTW <- list()
  for (i in 1:nbclusterDTW) {
    t <- summaryDTW[[i]]
    t <- gather(t, key = variable, value = value,-index)
    p <-
      ggplot(data = t, aes(x = index, y = value, color = variable))
    p <-
      p + geom_line() + labs(x = "Jour",
                             y = "Valeur",
                             title = paste("DTW, PAM, cluster", i))
    
    # Affichage du graphique avec plotly
    plotPAMDTW[[i]] <- ggplotly(p)
  }
  plotPAMSDTW <- list()
  for (i in 1:nbclusterSDTW) {
    t <- summarySDTW[[i]]
    t <- gather(t, key = variable, value = value,-index)
    p <-
      ggplot(data = t, aes(x = index, y = value, color = variable))
    p <-
      p + geom_line() + labs(
        x = "Jour",
        y = "Valeur",
        title = paste("Soft-DTW, PAM, cluster", i)
      )
    
    # Affichage du graphique avec plotly
    plotPAMSDTW[[i]] <- ggplotly(p)
  }
  
  plot2 <- subplot(plotPAMDTW, nrows = 2)
  plot3 <- subplot(plotPAMSDTW, nrows = 2)
  
})
tp <- tp["elapsed"]
print(paste("Les temps pour le calcul des matrice est de ", tp / 60, "minutes"))

print(plot1)
print(plot2)
print(plot3)
