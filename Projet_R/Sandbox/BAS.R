library(imputeTS)
library(readr)
library(ggplot2)
library(dtw)
library(DTWBI)
library(gridExtra)
library(kableExtra)
library(Rcpp)
source("EC_completion.r")

# P <- c(1, 4, 5, 10, 9, 3)
# P <- c(1, 7, 3, 4, 1, 10, 5)
# P <- c(2, 1, 4, 5, 1, 7, 3, 4, 1, 10, 5)
# P <- c(3, 1, 4, 5, 10, 9)

# Q <- c(1, 7, 3, 4, 1, 10, 5)

data <- tsNH4Complete

temp <- gapCreation(data,1000/length(data))$output_vector

P <- c(1,10,16,5)
Q <- c(1,4,13,9)

test <- dtw(x = P, y = Q, keep = TRUE)

print(test$costMatrix)

dtwPlot(test)

dtwPlotTwoWay(test)

dtwPlotThreeWay(test)

dtwPlotDensity(test)

# Prochains plot fait à partir de la demo de la librairie

# Asymetric step
ita <- dtw(P,Q,keep=TRUE,step=asymmetric)
plot(ita,type="density",main="Sine and cosine, asymmetric step")

## Windowing functions (global constraints) can be applied and plotted
dtwWindow.plot(sakoeChibaWindow, window.size = 10, main="So-called Itakura parallelogram window")

## Symmetric step with global parallelogram-shaped constraint
## Note how long (>2 steps) horizontal stretches are allowed within the window.
dtw(P,Q,keep=TRUE,window=itakuraWindow)->ita
dtwPlot(ita,type="density",main="Symmetric step with Itakura parallelogram window")

## Asymmetric step with slope constraint
## The "Itakura parallelogram" arises from the local constraint
## plus the boundary condition. Three sides of the parallelogram are seen
dtw(P,Q,keep=TRUE,step=typeIIIc)->ita
dtwPlot(ita,type="density",main="Slope-limited asymmetric step (Itakura)")

## Local and global constraint can be in effect at the same time
## Sakoe-chiba band, plus asymmetric step pattern
asyband<-dtw(P,Q,keep=TRUE,
             step=asymmetric,
             window.type=sakoeChibaWindow,
             window.size=30                  )

dtwPlot(asyband,type="density",main="Sine/cosine: asymmetric step, S-C window")