# Library
library(reticulate)

rm(list = ls())
#-------------------------------------------------------------------------------
# TimeGan
dataset <- datasets::co2
dataset <- r_to_py(dataset)

source_python("Python_function/TimeWarp.py")

rep <-
  TimeWarp(
    data = dataset,
    nbpts = 300,
    seq_len = 300,
    condition = 7,
    verbose = TRUE
  )

rep <- unlist(rep)

plot(rep,type = "l")

#-------------------------------------------------------------------------------
# Finding Feasible Windows


#-------------------------------------------------------------------------------
# Calculation of DTW and Soft-DTW distance matrices


#-------------------------------------------------------------------------------
# Clustering according to PAM