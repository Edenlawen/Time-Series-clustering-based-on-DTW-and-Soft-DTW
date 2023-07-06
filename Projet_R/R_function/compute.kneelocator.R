library(stats)

knee_locator <- function(x, y, direction = c("increasing", "decreasing", "both")) {
  n <- length(x)
  dx <- diff(x)
  dy <- diff(y)
  slopes <- dy / dx
  
  if (direction == "increasing") {
    slopes <- -slopes
  } else if (direction == "both") {
    slopes <- abs(slopes)
  }
  
  curve_lengths <- sqrt(dx^2 + dy^2)
  cumulative_lengths <- c(0, cumsum(curve_lengths))
  total_length <- cumulative_lengths[length(cumulative_lengths)]
  distances <- cumulative_lengths / total_length
  
  if (direction == "increasing") {
    distances <- 1 - distances
  }
  
  normalized_distances <- distances / max(distances)
  normalized_curvatures <- sqrt((y - predict(loess(y ~ x)))^2 + (x - x)^2) / max(distances)
  
  curvature_values <- normalized_distances + normalized_curvatures
  knee_index <- which.min(curvature_values)
  knee_x <- x[knee_index]
  knee_y <- y[knee_index]
  
  result <- list(knee_x = knee_x, knee_y = knee_y, knee_index = knee_index)
  class(result) <- "knee_locator_result"
  return(result)
}

get_knee <- function(knee_locator_result) {
  return(knee_locator_result$knee_x)
}
