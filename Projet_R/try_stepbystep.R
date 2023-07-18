# Library
library(reticulate)
library(readr)

rm(list = ls())
#-------------------------------------------------------------------------------
# TimeGan
dataset <- datasets::co2
dataset <- r_to_py(dataset)

source_python("Python_function/TimeWarp.py")

rep <-
  TimeWarp(
    data = dataset,
    nbpts = 50000,
    seq_len = 467,
    condition = 8,
    verbose = FALSE
  )

rep <- unlist(rep)

plot(rep,type = "l")

#-------------------------------------------------------------------------------
# Finding Feasible Windows
for (i in 1:25) {
  print(sample(1:10,1))
}


#-------------------------------------------------------------------------------
# Calculation of DTW and Soft-DTW distance matrices


#-------------------------------------------------------------------------------
# Clustering according to PAM