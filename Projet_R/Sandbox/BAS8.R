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
source("R_function/EC_completion.r")
source("R_function/EC_compareCourbe.r")
source("R_function/deformation.r")
source("R_function/soft_DTW_symmetric2.r")
source("R_function/globalF.r")


#' Algo de BAS8:
#'
#' Charger dataset @google
#' Initialiser une variable de nom @gapTaille ,
#' Initialiser une variable de nom @gapStart ,
#' Créer une copie du dataset contenant un gap de taille @gapTaille à un endroit défini par @gapStart
#' Initialiser une variable de nom @fenetresViable , ce sera un dataframe. Plus tard, on lui ajoutera les fenêtre acceptable d'après notre critère
#' Initialiser une variable de nom @queryTaille , elle permettra de savoir la taille de la query qu'on utilisera pour comparer
#' Initialiser une variable de nom @queryRef , la variable sera un vector et ayant le résultat de la fonction "globalfeatures" sur la query de reference
#' Initialiser une variable de nom @debut , elle sera égal à 1
#' Initialiser une variable de nom @fin , elle sera égal à @début + @queryTaille
#'
#' Tant que (@fin n'est pas égal à la taille du dataset)
#'    Si (!(@debut >= @gapStart ET @debut <= @gapStart + @gapTaille))
#'        Initialiser une variable de type vector et de nom @queryTemp , la variable sera un vector et ayant le résultat de la fonction "globalfeatures" pour  les données se trouvant dans le dataset google entre la position [@début;@fin]
#'        Initialiser une variable de nom @cosCompare , elle aura le résultat de la valeur absolue de la fonction "cosine" entre le vecteur @queryRef et @queryTemp telle que: abs(cosine(@queryRef , @queryTemp))
#'        Si @cosCompare >= 0.95
#'            On va ajouter @queryTemp à la variable @fenetresViable
#'        Fin Si
#'    Fin Si
#'    @debut += 1
#'    @fin += 1
#' Fin Tant que
#'
#' Initialiser une variable de type matrice et de nom @matriceDTW, elle sera de taille du: nombre de fenêtres retenus * nombre de fenêtres retenus
#' Boucle for @i de 1 au nombre de fenêtres retenus
#'    Boucle for @j de i jusqu'au nombre de fenêtres retenus
#'        @matriceDTW[i,j] sera égal au cout final de l'algo DTW entre les données à la position [,i] et la position [,j] de la variable @fenetreViable
#'        @matriceDTW[j,i] sera égal au cout final de l'algo DTW entre les données à la position [,i] et la position [,j] de la variable @fenetreViable
#'    Fin Boucle for
#' Fin Boucle for
#'
#' Regarder guideline.txt POUR LA SUITE, HISTOIRE DE KMEANS
#'

#-------------------------------------------------------------------------------
# Partie Code
tpsT <- system.time({

# data("google")
# data("co2")
# google <- co2
generated_data <- read_csv("generated_data.csv")
google <- as.ts(generated_data$`3.444022592320329750e+02`)

# Partie Chercher les fenêtres viables pour une taille de query donnée
tps <- system.time({
  gapTaille <- 5
  gapStart <- 400
  googleModif <-
    gapCreation(google, gapTaille / length(google), gapStart)$output_vector
  # queryTaille <- 10
  queryTaille <- 100
  featureRef <-
    globalfeatures(googleModif[(gapStart - queryTaille):(gapStart - 1)])
  queryRef <- googleModif[(gapStart - queryTaille):(gapStart - 1)]
  fenetresViable <- data.frame("queryRef" = queryRef)
  debut <- 1
  fin <- debut + queryTaille
  
  
  while (fin != length(google)) {
    if (!(debut %in% seq(gapStart - queryTaille, gapStart + gapTaille))) {
      featureTemp <- globalfeatures(googleModif[debut:(fin - 1)])
      queryTemp <- googleModif[debut:(fin - 1)]
      cosCompare <- abs(cosine(featureRef, featureTemp))
      
      if (!is.na(cosCompare) && cosCompare >= 0.995) {
        # print(debut)
        # print(cosCompare)
        fenetresViable <- cbind(fenetresViable,
                                queryTemp)
        colnames(fenetresViable)[ncol(fenetresViable)] <-
          paste0("Debut = ", debut)
      }
    }
    debut <- debut + 1
    fin <- fin + 1
  }
  rm(debut,fin,cosCompare)
})
tps <- tps["elapsed"]
print(paste(
  "Les temps pour la recherche de fenêtre viable est de ",
  tps,
  "secondes"
))

# Partie Calcul de la matrice DTW et Soft-DTW entre toutes les fenêtres viables
print("Calcul des matrices")
tps <- system.time({
  g = 0.001
  
  matriceDTW <-
    matrix(nrow = length(fenetresViable),
           ncol = length(fenetresViable))
  matriceSDTW <-
    matrix(nrow = length(fenetresViable),
           ncol = length(fenetresViable))
  print(length(fenetresViable))
  for (i in 1:length(fenetresViable)) {
    if(i%%50 == 0){print(paste("i=",i))}
    for (j in i:length(fenetresViable)) {
      matriceDTW[i, j] <-
        dtw_basic(fenetresViable[, i], fenetresViable[, j])
      matriceDTW[j, i] <- matriceDTW[i, j]
      matriceSDTW[i, j] <-
        sdtw(fenetresViable[, i], fenetresViable[, j], g)
      matriceSDTW[j, i] <- matriceSDTW[i, j]
    }
  }
})
tps <- tps["elapsed"]
print(paste("Les temps pour le calcul des matrice est de ", tps, "secondes"))


# Partie kmeans
print("Début de la partie Kmeans")
nbclusterDTW <-
  fviz_nbclust(matriceDTW, kmeans, method = "silhouette")
nbclusterDTW <- nbclusterDTW$data$y
nbclusterDTW <- which.max(nbclusterDTW)
nbclusterSDTW <-
  fviz_nbclust(matriceSDTW, pam, method = "silhouette")
nbclusterSDTW <- nbclusterSDTW$data$y
nbclusterSDTW <- which.max(nbclusterSDTW)

matriceDTWdist <- dist(matriceDTW)
matriceSDTWdist <- dist(matriceSDTW)

resultatKmeansDTW <- kmeans(matriceDTWdist, nbclusterDTW)
resultatKmeansSDTW <- kmeans(matriceSDTWdist, nbclusterSDTW)

plot1 <-
  fviz_cluster(resultatKmeansDTW, data = matriceDTWdist) + ggtitle("Kmeans avec DTW")
plot2 <-
  fviz_cluster(resultatKmeansSDTW, data = matriceSDTWdist) + ggtitle("Kmeans avec Soft-DTW")

#   Partie 1: Cluster ayant la moyenne DTW/SDTW la plus faible
avgClusterDTW <- rep(0, times = nbclusterDTW)
avgClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 2:length(fenetresViable)) {
  avgClusterDTW[resultatKmeansDTW$cluster[i]] <-
    avgClusterDTW[resultatKmeansDTW$cluster[i]] + matriceDTWdist[i]
  avgClusterSDTW[resultatKmeansSDTW$cluster[i]] <-
    avgClusterSDTW[resultatKmeansSDTW$cluster[i]] + matriceSDTWdist[i]
}

for (i in 1:nbclusterSDTW) {
  avgClusterDTW[i] <-
    avgClusterDTW[i] / table(resultatKmeansDTW$cluster)[i]
  avgClusterSDTW[i] <-
    avgClusterSDTW[i] / table(resultatKmeansSDTW$cluster)[i]
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
  resultatKmeansDTW$cluster
)), "a le plus de points"))
print(paste("Pour SDTW, le cluster", which.max(table(
  resultatKmeansSDTW$cluster
)), "a le plus de points"))


# Afficher graph par cluster
dfkmeanDTW <- list()
for (i in 1:nbclusterDTW) {
  dfkmeanDTW <-
    append(dfkmeanDTW, list(data.frame("index" = 1:queryTaille)), after = length(dfkmeanDTW))
}
dfkmeanSDTW <- list()
for (i in 1:nbclusterSDTW) {
  dfkmeanSDTW <-
    append(dfkmeanSDTW, list(data.frame("index" = 1:queryTaille)), after = length(dfkmeanSDTW))
}
for (i in 2:length(fenetresViable)) {
  df_selected <- dfkmeanDTW[[resultatKmeansDTW$cluster[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfkmeanDTW[[resultatKmeansDTW$cluster[i]]] <- df_selected
  df_selected <- dfkmeanSDTW[[resultatKmeansSDTW$cluster[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfkmeanSDTW[[resultatKmeansSDTW$cluster[i]]] <- df_selected
}
plotkmeansDTW <- list()
for (i in 1:nbclusterDTW) {
  t <- dfkmeanDTW[[i]]
  t <- gather(t, key = variable, value = value,-index)
  p <-
    ggplot(data = t, aes(x = index, y = value, color = variable))
  p <-
    p + geom_line() + labs(x = "Jour",
                           y = "Valeur",
                           title = paste("DTW, Kmeans, cluster", i))
  
  # Affichage du graphique avec plotly
  plotkmeansDTW[[i]] <- ggplotly(p)
}
plotkmeansSDTW <- list()
for (i in 1:nbclusterDTW) {
  t <- dfkmeanSDTW[[i]]
  t <- gather(t, key = variable, value = value,-index)
  p <-
    ggplot(data = t, aes(x = index, y = value, color = variable))
  p <-
    p + geom_line() + labs(x = "Jour",
                           y = "Valeur",
                           title = paste("Soft-DTW, Kmeans, cluster", i))
  
  # Affichage du graphique avec plotly
  plotkmeansSDTW[[i]] <- ggplotly(p)
}

print(plotkmeansDTW)
print(plotkmeansSDTW)

print("Fin de la partie Kmeans")

#-------------------------------------------------------------------------------
# Partie PAM (K-medois)
print("Début de la partie PAM")

nbclusterDTW <-
  fviz_nbclust(matriceDTW, kmeans, method = "silhouette")
nbclusterDTW <- nbclusterDTW$data$y
nbclusterDTW <- which.max(nbclusterDTW)
nbclusterSDTW <-
  fviz_nbclust(matriceSDTW, pam, method = "silhouette")
nbclusterSDTW <- nbclusterSDTW$data$y
nbclusterSDTW <- which.max(nbclusterSDTW)

resultatPamDTW <- pam(matriceDTWdist, nbclusterDTW)
resultatPamSDTW <- pam(matriceSDTWdist, nbclusterSDTW)

# Réduction de dimension avec PCA
pcaDTW <- prcomp(matriceDTWdist)
pcaSDTW <- prcomp(matriceSDTWdist)

# Utilisation des coordonnées PCA pour visualiser les clusters
plot3 <-
  fviz_pca_ind(pcaDTW, habillage = as.factor(resultatPamDTW$clustering)) + ggtitle("K-Médoïdes avec DTW")
plot4 <-
  fviz_pca_ind(pcaSDTW, habillage = as.factor(resultatPamSDTW$clustering)) + ggtitle("K-Médoïdes avec Soft-DTW")

rm(pcaDTW,pcaSDTW)

#   Partie 1: Cluster ayant la moyenne DTW/SDTW la plus faible
avgClusterDTW <- rep(0, times = nbclusterDTW)
avgClusterSDTW <- rep(0, times = nbclusterSDTW)

for (i in 2:length(fenetresViable)) {
  avgClusterDTW[resultatPamDTW$clustering[i]] <-
    avgClusterDTW[resultatPamDTW$clustering[i]] + matriceDTWdist[i]
  avgClusterSDTW[resultatPamSDTW$clustering[i]] <-
    avgClusterSDTW[resultatPamSDTW$clustering[i]] + matriceSDTWdist[i]
}

for (i in 1:nbclusterSDTW) {
  avgClusterDTW[i] <-
    avgClusterDTW[i] / table(resultatPamDTW$clustering)[i]
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

# Afficher graph par cluster
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
for (i in 2:length(fenetresViable)) {
  df_selected <- dfPAMDTW[[resultatPamDTW$clustering[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfPAMDTW[[resultatPamDTW$clustering[i]]] <- df_selected
  df_selected <- dfPAMSDTW[[resultatPamSDTW$clustering[i]]]
  df_selected <- cbind(df_selected, fenetresViable[i])
  dfPAMSDTW[[resultatPamSDTW$clustering[i]]] <- df_selected
}
plotPAMDTW <- list()
for (i in 1:nbclusterDTW) {
  t <- dfPAMDTW[[i]]
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
for (i in 1:nbclusterDTW) {
  t <- dfPAMSDTW[[i]]
  t <- gather(t, key = variable, value = value,-index)
  p <-
    ggplot(data = t, aes(x = index, y = value, color = variable))
  p <-
    p + geom_line() + labs(x = "Jour",
                           y = "Valeur",
                           title = paste("Soft-DTW, PAM, cluster", i))
  
  # Affichage du graphique avec plotly
  plotPAMSDTW[[i]] <- ggplotly(p)
}

print(plotPAMDTW)
print(plotPAMSDTW)

print("Fin de la partie PAM")

# Partie affichage graphique
df <- data.frame(index = 1:length(google))
df$ref <- google
df$DTWBI <- compute.DTWBI_QF(googleModif, 6, 0, queryTaille)

ggplot(df, aes(x = index)) +
  geom_line(aes(y = ref, color = "Base")) +
  geom_line(aes(y = DTWBI, color = "SDTW")) +
  coord_cartesian(xlim = c(gapStart - 5, gapStart + gapTaille + 5))



# Création du dataframe de données
dfTEMP <- data.frame(x = 1:queryTaille)
dfTEMP <- bind_cols(dfTEMP, fenetresViable)

# Fusion des données avec tidyr
dfTEMP <- gather(dfTEMP, key = variable, value = value,-x)

# Création de l'objet ggplot
p <-
  ggplot(data = dfTEMP, aes(x = x, y = value, color = variable))

# Ajout des courbes
p <- p + geom_line() + labs(x = "Jour", y = "Valeur")

# Affichage du graphique avec plotly
plot <- ggplotly(p)

# Sauvegarde du graphique interactif en HTML
htmlwidgets::saveWidget(plot, "test.html")

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

})
tpsT <- tpsT["elapsed"]
print(paste("Les temps pour la recherche de fenêtre viable est de ",tpsT,"secondes"))
