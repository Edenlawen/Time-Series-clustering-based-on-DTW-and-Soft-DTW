library(TSA)
library(imputeTS)
library(DTWBI)
library(dtw)
library(ggplot2)
library(imputeTS)
library(readr)
library(gridExtra)
library(kableExtra)
library(Rcpp)
library(cluster)
library(factoextra)
library(tidyr)
library(dplyr)
library(kneedle)
source("R_function/EC_completion.r")
source("R_function/EC_compareCourbe.r")
source("R_function/deformation.r")

data("google")

AY <- gapCreation(google, 20 / length(google), 160)$output_vector

cond <- 0

score <- NULL
querys <- NULL

for (i in 4:40) {
  sig <-
    compute.DTWBI_QF(
      AY,
      acceptedHole = 40,
      smallHole = 9,
      S_Query = i,
      verbose = FALSE
    )
  valeur <- compute.similariteAire(google, sig)
  if (i %% 20 == 0) {
    print(i)
  }
  score <- c(score, valeur)
  if (cond < valeur) {
    query = i
    cond = valeur
    print(i)
    print("TOP")
    print(valeur)
    
  }
  if (valeur > 0.999975) {
    querys <- c(querys, i)
  }
}
res <- list(var1 = query, var2 = score)

df <- data.frame(compute.DTWBI_QF(
  AY,
  acceptedHole = 50,
  smallHole = 9,
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
      smallHole = 9,
      S_Query = querys[i],
      verbose = FALSE
    )
  )
  colnames(df)[ncol(df)] <- paste0("Query = ", querys[i])
}

dft <- t(df)
dft1 <- as.data.frame(dft)
dft2 <- as.data.frame(dft)

max_distinct <- max(apply(dft1, 2, function(col) length(unique(col))))

inertie_kmeans <- sapply(2:max_distinct, function(k) {
  kmeans_res <- kmeans(dft1, centers = k)
  kmeans_res$tot.withinss
})

kn <- kneedle(2:max_distinct,inertie_kmeans, decreasing = TRUE, concave = TRUE)
plot(2:4,inertie_kmeans, type = "b")
abline(v=kn)


kmea <- kmeans(dft1,
               centers = kn[1],
               iter.max = 20,
               nstart = 1)
dft1$cluster <- kmea$cluster
df_plot1 <- dft1 %>%
  gather(key = "Query", value = "Value", -cluster)
ggplot(df_plot1, aes(
  x = Query,
  y = Value,
  color = factor(cluster)
)) +
  geom_line() +
  geom_point() +
  labs(title = "Résultats du clustering k-means") +
  theme_minimal()

kmed <- pam(dft2, k = kn[1], metric = "euclidean")
dft2$cluster <- kmed$clustering
df_plot2 <- dft2 %>%
  gather(key = "Query", value = "Value", -cluster)
ggplot(df_plot2, aes(
  x = Query,
  y = Value,
  color = factor(cluster)
)) +
  geom_point() +
  labs(title = "Résultats du clustering PAM") +
  theme_minimal()

#-------------------------------------------------------------------------------
temp <- NULL
for (i in 1:length(google)) {
  sumtemp <- 0
  for (j in 1:length(df[i,])) {
    sumtemp <- sumtemp + df[i,j]
  }
  sumtemp <- sumtemp / length(df[i,])
  temp <- c(temp,sumtemp)
}
dfmean <- data.frame(temp, google)

ggplot(dfmean, aes(x = 1:521)) +
  geom_line(aes(y = temp, color = "DTW mean")) +
  geom_line(aes(y = google, color = "google")) +
  coord_cartesian(xlim = c(150, 180))

print(compute.indicateurComp(dfmean$google,dfmean$temp))

#-------------------------------------------------------------------------------
dfSim <- data.frame(res$var2)

nb_points_distincts <- dfSim %>%
  distinct(res.var2) %>%
  n_distinct()

inertie_kmeans <- sapply(2:nb_points_distincts, function(k) {
  kmeans_res <- kmeans(dfSim, centers = k)
  kmeans_res$tot.withinss
})


kn <- kneedle(2:nb_points_distincts,inertie_kmeans, concave = FALSE, decreasing = TRUE)
plot(2:nb_points_distincts,inertie_kmeans, type = "b")
abline(v=kn)


kmea <- kmeans(dfSim,
               centers = kn[1],
               iter.max = 20,
               nstart = 1)
dfSim$cluster <- kmea$cluster
df_plot3 <- dfSim %>%
  gather(key = "Query", value = "Value", -cluster)
ggplot(df_plot3, aes(
  x = Query,
  y = Value,
  color = factor(cluster)
)) +
  geom_point() +
  labs(title = "Résultats du clustering k-means") +
  theme_minimal()
