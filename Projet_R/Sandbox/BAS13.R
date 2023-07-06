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

source("R_function/EC_completion.r")
source("R_function/globalF.r")
source("R_function/compute.erreurRMSE.r")
source("R_function/EC_compareCourbe.r")

dtwresult <- rep(0, times = 3)

tp <- system.time({
  # google <- read_csv("df_filled_1W_730.csv")
  # google <- google$Fluo_FFU[1:20000]
  
  google <- datasets::co2
  
  gan <- read_csv("generated_data.csv")
  google <- c(google,gan$`3.444022592320329750e+02`)
  
  g = 0.001
  gapTaille <- 7
  compteur <- 1
  queryTaille <- 12
  boucle <- length(datasets::co2) - 25
  
  # for (m in sample(50:19950, 20000, replace = TRUE)[6000:(6000+boucle)]) {
  for (m in 25:boucle) {
    print(paste(((compteur*100)/boucle),"%"))
    gapStart <- m
    dataModif <-
      gapCreation(google, gapTaille / length(google), gapStart)$output_vector
    featureRef <-
      globalfeatures(dataModif[(gapStart - queryTaille):(gapStart - 1)])
    queryRef <- dataModif[(gapStart - queryTaille):(gapStart - 1)]
    fenetresViable <- data.frame("queryRef" = queryRef)
    repRef <- google[gapStart:(gapStart + gapTaille - 1)]
    reponseVialbe <- data.frame("repRef" = repRef)
    debut <- 1
    fin <- debut + queryTaille
    
    if (length(google) < 1000) {
      step_threshold = 2
    } else{
      if (length(google) > 10000) {
        step_threshold = 50
      } else{
        step_threshold = 10
      }
    }
    
    if (length(google) < 10000) {
      threshold_cos = 0.9995
    } else {
      threshold_cos = 0.995
    }
    
    
    stepBase <- step_threshold
    while (fin <= (length(google) - queryTaille - gapTaille)) {
      if (!(
        debut %in% seq(
          gapStart - queryTaille - gapTaille,
          gapStart + queryTaille + gapTaille
        )
      )) {
        featureTemp <- globalfeatures(dataModif[debut:(fin - 1)])
        queryTemp <- dataModif[debut:(fin - 1)]
        cosCompare <- abs(cosine(featureRef, featureTemp))
        
        if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
          fenetresViable <- cbind(fenetresViable,
                                  queryTemp)
          colnames(fenetresViable)[ncol(fenetresViable)] <-
            paste0("Debut = ", debut)
          
          repTemp <- dataModif[fin:(fin + gapTaille - 1)]
          
          reponseVialbe <- cbind(reponseVialbe, repTemp)
          colnames(reponseVialbe)[ncol(reponseVialbe)] <-
            paste0("Debut = ", debut)
        }
      }
      debut <- debut + step_threshold
      fin <- fin + step_threshold
      
      if (fin >= (length(google) - queryTaille - gapTaille)) {
        if (length(fenetresViable) < 200) {
          if (step_threshold == 1) {
            threshold_cos <- 0.95
            step_threshold <- stepBase
          }
          debut <- 1
          fin <- debut + queryTaille
          step_threshold <-
            step_threshold - floor(step_threshold / 2)
          # print(step_threshold)
          rm(fenetresViable, reponseVialbe)
          fenetresViable <- data.frame("queryRef" = queryRef)
          reponseVialbe <- data.frame("repRef" = repRef)
        }
      }
    }
    fenetresViable <- subset(fenetresViable, select = -1)
    reponseVialbe <- subset(reponseVialbe, select = -1)
    rm(debut, fin, cosCompare, featureRef, queryTemp)
    
    # print("fin fenetre")
    
    # Partie Calcul de la matrice DTW et Soft-DTW entre toutes les fenêtres viables
    
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
    
    if (length(fenetresViable) < 25) {
      kmax <- length(fenetresViable) - 1
    } else{
      kmax <- 25
    }
    
    nbclusterDTW <-
      fviz_nbclust(
        t(fenetresViable),
        pam,
        diss = matriceDTW,
        method = "silhouette",
        k.max = kmax
      )
    nbclusterDTW <- nbclusterDTW$data$y
    nbclusterDTW <- which.max(nbclusterDTW)
    nbclusterSDTW <-
      fviz_nbclust(
        t(fenetresViable),
        pam,
        diss = matriceSDTW,
        method = "silhouette",
        k.max = kmax
      )
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
    
    # cat("\nPartie 1\n")
    # print(paste(
    #   "Pour DTW, le cluster",
    #   which.min(avgClusterDTW),
    #   "a la moyenne de coût DTW la plus faible"
    # ))
    # print(paste(
    #   "Pour SDTW, le cluster",
    #   which.min(avgClusterSDTW),
    #   "a la moyenne de coût DTW la plus faible"
    # ))
    
    # Partie 5:
    
    avgAmpClusterDTW <- rep(0, times = nbclusterDTW)
    avgAmpClusterSDTW <- rep(0, times = nbclusterSDTW)
    
    for (i in 1:length(reponseVialbe)) {
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
    
    # cat("\nPartie 5\n")
    # print(
    #   paste(
    #     "Pour DTW, le cluster",
    #     which.min(avgAmpClusterDTW),
    #     "a la moyenne de coût DTW des réponses la plus faible"
    #   )
    # )
    # print(
    #   paste(
    #     "Pour SDTW, le cluster",
    #     which.min(avgAmpClusterSDTW),
    #     "a la moyenne de coût DTW des réponses la plus faible"
    #   )
    # )
    
    repC1DTW <-
      data.frame("repRef" = google[gapStart:(gapStart + gapTaille - 1)])
    repC1SDTW <-
      data.frame("repRef" = google[gapStart:(gapStart + gapTaille - 1)])
    
    repC5DTW <-
      data.frame("repRef" = google[gapStart:(gapStart + gapTaille - 1)])
    repC5SDTW <-
      data.frame("repRef" = google[gapStart:(gapStart + gapTaille - 1)])
    
    for (i in 1:length(reponseVialbe)) {
      if (resultatPamDTW[i] == which.min(avgClusterDTW)) {
        repC1DTW <- cbind(repC1DTW, reponseVialbe[, i])
      }
      if (resultatPamSDTW[i] == which.min(avgClusterSDTW)) {
        repC1SDTW <- cbind(repC1SDTW, reponseVialbe[, i])
      }
      if (resultatPamDTW[i] == which.min(avgAmpClusterDTW)) {
        repC5DTW <- cbind(repC5DTW, reponseVialbe[, i])
      }
      if (resultatPamSDTW[i] == which.min(avgAmpClusterSDTW)) {
        repC5SDTW <- cbind(repC5SDTW, reponseVialbe[, i])
      }
    }
    repC1DTW <- subset(repC1DTW, select = -1)
    repC1SDTW <- subset(repC1SDTW, select = -1)
    repC5DTW <- subset(repC5DTW, select = -1)
    repC5SDTW <- subset(repC5SDTW, select = -1)
    
    repC1DTW <- t(repC1DTW)
    repC1SDTW <- t(repC1SDTW)
    repC5DTW <- t(repC5DTW)
    repC5SDTW <- t(repC5SDTW)
    
    
    meanRepC1DTW <- vector("numeric", length = 0)
    meanRepC1SDTW <- vector("numeric", length = 0)
    
    meanRepC5DTW <- vector("numeric", length = 0)
    meanRepC5SDTW <- vector("numeric", length = 0)
    
    for (i in 1:gapTaille) {
      # meanRepC1DTW <- c(meanRepC1DTW, mean(repC1DTW[, i]))
      # meanRepC1SDTW <- c(meanRepC1SDTW, mean(repC1SDTW[, i]))
      # 
      # meanRepC5DTW <- c(meanRepC5DTW, mean(repC5DTW[, i]))
      # meanRepC5SDTW <- c(meanRepC5SDTW, mean(repC5SDTW[, i]))
      
      meanRepC1DTW <- c(meanRepC1DTW, quantile(repC1DTW[, i], 0.5))
      meanRepC1SDTW <-
        c(meanRepC1SDTW, quantile(repC1SDTW[, i], 0.5))

      meanRepC5DTW <- c(meanRepC5DTW, quantile(repC5DTW[, i], 0.5))
      meanRepC5SDTW <-
        c(meanRepC5SDTW, quantile(repC5SDTW[, i], 0.5))
      
    }
    
    df <-
      data.frame(
        "index" = 1:length(meanRepC1DTW),
        "main" = google[gapStart:(gapStart + gapTaille - 1)],
        "avgC1DTW" = meanRepC1DTW,
        "avgC1SDTW" = meanRepC1SDTW,
        "avgC5DTW" = meanRepC5DTW,
        "avgC5SDTW" = meanRepC5SDTW
      )
    
    distfin <-
      which.min(c(
        abs(google[gapStart + gapTaille] - df[length(df), 3]),
        abs(google[gapStart + gapTaille] - df[length(df), 4]),
        abs(google[gapStart + gapTaille] - df[length(df), 5]),
        abs(google[gapStart + gapTaille] - df[length(df), 6])
      ))
    
    dtwbi <- compute.DTWBI_QF(dataModif, 20, 0, queryTaille, verbose = FALSE)[(gapStart):(gapStart+gapTaille-1)]
    dataModif[gapStart:(gapStart + gapTaille - 1)] <- df[, distfin+2]
    
    df2 <- data.frame(
      "index" = 1:gapTaille,
      "main" = google[(gapStart):(gapStart+gapTaille-1)],
      "comp" = dataModif[(gapStart):(gapStart+gapTaille-1)],
      "DTWBI" = dtwbi
    )
    
    temp <- rep(0, times = 2)
    
    temp[1] <- compute.indicateurComp(df2$main,df2$comp)[12]
    temp[2] <- compute.indicateurComp(df2$main,df2$DTWBI)[12]
    
    if(temp[1]==temp[2]){
      dtwresult[3] <- dtwresult[3] + 1
      print(m)
    }else{
      dtwresult[which.max(temp)] <- dtwresult[which.max(temp)] + 1
    }
    rm(dataModif, fenetresViable, reponseVialbe)
    compteur <- compteur + 1
  }
})
tp <- tp["elapsed"]
print(paste("Les temps pour le programme est de ", tp / 60, "minutes"))

print(dtwresult)
cat("\nMa méthode:\n")
print(dtwresult[1])
cat("\nDTWBI:\n")
print(dtwresult[2])
cat("\nSame:\n")
print(dtwresult[3])