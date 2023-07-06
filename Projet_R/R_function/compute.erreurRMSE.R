rmse <- function(signal1, signal2) {
  if (length(signal1) != length(signal2)) {
    stop("Les signaux doivent avoir la même longueur.")
  }
  
  squared_diff <- (signal1 - signal2)^2
  mean_squared_error <- mean(squared_diff)
  rmse_value <- sqrt(mean_squared_error)
  
  return(rmse_value)
}