# Library
library(reticulate)
library(readr)
library(ggplot2)
library(plotly)
library(tidyr)
library(reshape2)
library(tidyverse)

rm(list = ls())
#-------------------------------------------------------------------------------
# TimeGan
dataset <- datasets::co2
dataset <- r_to_py(dataset)

# dataset <- read_csv("csv/bochoiPhuLien.csv")
# dataset <- dataset[,-1]
# temp <- unlist(t(dataset))
# temp <- as.data.frame(temp)

# temp$index <- 1:12
# 
# date <- 1959
# for (i in 1:(2015 - 1959 + 1)) {
#   colnames(temp)[i] <- date
#   date <- date + 1
# }
# 
# t <- gather(temp, key = variable, value = value, -index)
# p <-
#   ggplot(data = t, aes(x = index, y = value, color = variable))
# p <-
#   p + geom_line() + labs(x = "Jour",
#                          y = "Valeur")
# 
# # Affichage du graphique avec plotly
# print(ggplotly(p))

# dataset <- as.vector(t(dataset))
# dataset <- na.omit(dataset)

# dataset <- r_to_py(dataset)

source_python("Python_function/TimeWarp.py")

rep <-
  TimeWarp(
    data = dataset,
    nbpts = 20000,
    seq_len = 467,
    condition = 8,
    verbose = FALSE
  )

rep <- unlist(rep)

plot(rep, type = "l")
abline(v = 467, col = "red", lwd = 3)

#-------------------------------------------------------------------------------
# Finding Feasible Windows
source("R_function/compute.searchWindows.r")

fenetresViable <- compute.searchWindows(rep, rep[400:407], 100, 250, random = TRUE, verbose = FALSE)


#-------------------------------------------------------------------------------
# Calculation of DTW and Soft-DTW distance matrices
source("R_function/compute.DistanceMatrix.r")

infoMatrixDTW <- compute.DistanceMatrixDTW(fenetresViable)
infoMatrixSDTW <- compute.DistanceMatrixSDTW(fenetresViable, 0.001)

#-------------------------------------------------------------------------------
# Clustering according to PAM
source("R_function/compute.PAMClustering.r")

PAMDTW <-
  compute.PAMClustering(fenetresViable, infoMatrixDTW[[1]], queryRef = rep[400:407], nbcluster = NULL)
PAMSDTW <-
  compute.PAMClustering(fenetresViable, infoMatrixSDTW[[1]], queryRef = rep[400:407], nbcluster = NULL)
