library(Rcpp)
library(ggplot2)
library(dtw)
library(DTWBI)
source("R_function/EC_completion.r")
source("R_function/fonction_maison.r")
source("R_function/EC_compareCourbe.r")

sourceCpp("Cpp_function/dtw.cpp")

# P <- c(1,4,13,9)
# Q <- c(1,10,16,5)

# P <- c(1,4,16,8,2,1)
# Q <- c(1,5,10,15,20,18)

P <- c(1, 4, 5, 10, 9, 3, 5)
Q <- c(1, 4, 5, 10, 8, 3, 5)
# Q <- c(1, 7, 3, 4, 1, 10, 5)

test <- costMatrix(P,Q)
print(test)

test2 <- dtw(x = P, y = Q, keep.internals = TRUE,  distance.only = FALSE, 
             step.pattern = symmetric2)
print(test2$costMatrix)

## Well-known step patterns
# symmetric1
# symmetric2
# asymmetric

## Slope-constrained step patterns from Sakoe-Chiba (Sakoe1978)
# symmetricP0;  asymmetricP0
# symmetricP05; asymmetricP05
# symmetricP1;  asymmetricP1
# symmetricP2;  asymmetricP2

## Step patterns classified according to Rabiner-Myers (Myers1980)
# typeIa;   typeIb;   typeIc;   typeId;
# typeIas;  typeIbs;  typeIcs;  typeIds;  # smoothed
# typeIIa;  typeIIb;  typeIIc;  typeIId;
# typeIIIc; typeIVc;

## Miscellaneous
# mori2006;
# rigid;