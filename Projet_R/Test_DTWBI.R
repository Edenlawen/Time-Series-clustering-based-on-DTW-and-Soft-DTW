library(DTWBI)
library(readr)
library(imputeTS)

data("google")
query <- google
query <- as.numeric(query)
ref <- google
ref <- as.numeric(ref)

query <- gapCreation(query, rate = 0.1)

t_gap <- query$begin_gap
T_gap <- query$gap_size
results_DTWBI <- DTWBI_univariate(query$output_vector, t_gap, T_gap)


plot(results_DTWBI$output_vector, type = "l", col = "blue", pch = 16)
#lines(results_DTWBI$output_vector, type="o", col="red")

######################################################################
par(mfrow = c(1, 2))

plot(results_DTWBI$input_vector, type = "o", col = "blue", pch = 16)
lines(results_DTWBI$input_vector, type="o", col="red")

plot(results_DTWBI$output_vector, type = "o", col = "blue", pch = 16)
lines(results_DTWBI$output_vector, type="o", col="red")

######################################################################
inter1 <- na_interpolation(query$output_vector, option = "linear")
inter2 <- na_interpolation(query$output_vector, option = "spline")
inter3 <- na_interpolation(query$output_vector, option = "stine")

a <- 270
b <- 330

par(mfrow = c(2, 2))

plot(results_DTWBI$output_vector[a:b], type = "o", col = "blue",
     pch = 16, main = "DTWBI")
lines(results_DTWBI$output_vector[a:b], type="o", col="red")

plot(inter1[a:b], type = "o", col = "blue",
     pch = 16, main = "na_interpolation de type linéaire")
lines(inter1[a:b], type="o", col="red")

plot(inter2[a:b], type = "o", col = "blue",
     pch = 16, main = "na_interpolation de type Spline")
lines(inter2[a:b], type="o", col="red")

plot(inter3[a:b], type = "o", col = "blue",
     pch = 16, main = "na_interpolation de type Stine")
lines(inter3[a:b], type="o", col="red")