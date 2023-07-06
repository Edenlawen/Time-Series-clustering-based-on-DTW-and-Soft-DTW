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
library(Rcpp)
source("R_function/EC_completion.r")
# source("R_function/EC_compareCourbe.r")
source("R_function/globalF.r")
sourceCpp("Cpp_function/dtw.cpp")

google <- read_csv("generated_data.csv")
google <- google$`3.444022592320329750e+02`

# Partie Chercher les fenêtres viables pour une taille de query donnée
print("Chercher fenetre viable")
tps <- system.time({
  gapTaille <- 5
  gapStart <- 400
  dataModif <-
    gapCreation(google, gapTaille / length(google), gapStart)$output_vector
  queryTaille <- 12
  featureRef <-
    globalfeatures(dataModif[(gapStart - queryTaille):(gapStart - 1)])
  queryRef <- dataModif[(gapStart - queryTaille):(gapStart - 1)]
  fenetresViable <- data.frame("queryRef" = queryRef)
  debut <- 1
  fin <- debut + queryTaille
  
  if (length(google) < 1000) {
    step_threshold = 2
  }
  else{
    if (length(google) > 10000) {
      step_threshold = 50
    }
    else{
      step_threshold = 10
    }
  }
  
  if (length(google) < 10000) {
    threshold_cos = 0.9995
  } else {
    threshold_cos = 0.995
  }
  
  
  while (fin <= length(google)) {
    if (!(debut %in% seq(gapStart - queryTaille, gapStart + gapTaille))) {
      featureTemp <- globalfeatures(dataModif[debut:(fin - 1)])
      queryTemp <- dataModif[debut:(fin - 1)]
      cosCompare <- abs(cosine(featureRef, featureTemp))
      
      #if (!is.na(cosCompare) && cosCompare >= 0.9999999) {
      if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
        # print(debut)
        # print(cosCompare)
        fenetresViable <- cbind(fenetresViable,
                                queryTemp)
        colnames(fenetresViable)[ncol(fenetresViable)] <-
          paste0("Debut = ", debut)
      }
    }
    debut <- debut + step_threshold
    fin <- fin + step_threshold
  }
  rm(debut, fin, cosCompare, featureRef, queryTemp)
})
tps <- tps["elapsed"]
print(paste(
  "Les temps pour la recherche de fenêtre viable est de ",
  tps,
  "secondes"
))

print("Calcul des matrices")
tps <- system.time({
  g = 0.001
  
  matriceDTW <-
    matrix(0,
           nrow = length(fenetresViable),
           ncol = length(fenetresViable))
  print(length(fenetresViable))
  for (i in 1:length(fenetresViable)) {
    if (i %% 100 == 0) {
      print(paste("i=", i))
    }
    for (j in i:length(fenetresViable)) {
      # matriceDTW[j, i] <- costFin(fenetresViable[, i], fenetresViable[, j])
      matriceDTW[j, i] <-
        dtw_basic(fenetresViable[, i], fenetresViable[, j])
    }
  }
  matriceDTW <-
    t(matriceDTW) + matriceDTW - diag(diag(matriceDTW))
  # write.csv(matriceDTW, file = "matriceDTW.csv", row.names = FALSE)
  # write.csv(matriceSDTW, file = "matriceSDTW.csv", row.names = FALSE)
})
tps <- tps["elapsed"]
print(tps)
rm(g, i, j)

# sig1 <- c(2,8,3)
# sig2 <- c(2,6,8)
#
# temp <- costFin(sig1,sig2)

# Autre facon de calculer la matrice de distance
# proxy::dist(t(fenetresViable), method = "dtw_basic", upper = FALSE, diag = TRUE)

bo <-
  proxy::dist(
    t(fenetresViable),
    method = "dtw",
    window.type = "itakura",
    upper = FALSE,
    diag = TRUE
  )
bo1 <-
  proxy::dist(t(fenetresViable),
              method = "dtw",
              upper = FALSE,
              diag = FALSE)
bo2 <-
  proxy::dist(t(fenetresViable),
              method = "dtw_basic",
              upper = FALSE,
              diag = TRUE)
bo3 <-
  proxy::dist(
    t(fenetresViable),
    method = "dtw",
    window.type = "sakoechiba",
    window.size = 5,
    upper = FALSE,
    diag = TRUE
  )
bo4 <-
  proxy::dist(
    t(fenetresViable),
    method = "sdtw",
    gamma = 0.001,
    upper = FALSE,
    diag = TRUE
  )

subplot(
  fviz_nbclust(t(fenetresViable), pam, method = "silhouette", k.max = 25),
  fviz_nbclust(
    t(fenetresViable),
    pam,
    diss = bo,
    method = "silhouette",
    k.max = 25
  ),
  fviz_nbclust(
    t(fenetresViable),
    pam,
    diss = bo1,
    method = "silhouette",
    k.max = 25
  ),
  fviz_nbclust(
    t(fenetresViable),
    pam,
    diss = bo2,
    method = "silhouette",
    k.max = 25
  ),
  fviz_nbclust(
    t(fenetresViable),
    pam,
    diss = bo3,
    method = "silhouette",
    k.max = 25
  ),
  fviz_nbclust(
    t(fenetresViable),
    pam,
    diss = bo4,
    method = "silhouette",
    k.max = 25
  ),
  nrows = 2
)

avg <- 0
for (k in 1:7) {
  print(k)
  tps <- system.time({
    g = 0.001
    
    matriceDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "dtw",
        window.type = "sakoechiba",
        window.size = 10,
        upper = FALSE,
        diag = TRUE
      )
    matriceSDTW <-
      proxy::dist(
        t(fenetresViable),
        method = "sdtw",
        gamma = 0.001,
        upper = FALSE,
        diag = TRUE
      )
  })
  tps <- tps["elapsed"]
  avg <- avg + tps
}
avg <- avg / 7