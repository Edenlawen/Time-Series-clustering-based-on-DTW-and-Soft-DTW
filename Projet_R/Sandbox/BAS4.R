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
source("R_function/EC_completion.r")
source("R_function/EC_compareCourbe.r")
source("R_function/fonction_maison.r")

df_filled_1W_730 <- read_csv("df_filled_1W_730.csv")

ggplot(df_filled_1W_730, aes(x = time)) +
  geom_line(aes(y = Fluo_FFU))

infoNA(df_filled_1W_730, fig = F)

test <-
  compute.DTWBI_QF(
    df_filled_1W_730$Fluo_FFU,
    acceptedHole = 4725,
    smallHole = 0,
    S_Query = 24,
    verbose = FALSE
  )

ggplot(df_filled_1W_730, aes(x = time)) +
  geom_line(aes(y = Fluo_FFU, colour = "Main")) +
  geom_line(aes(y = test, colour = "DTWBI")) +
  scale_x_continuous(limits = c(df_filled_1W_730$time[80331], df_filled_1W_730$time[80441]))
#-------------------------------------------------------------------------------
f2 <- function(Liste,
               taille_gap,
               start_gap,
               verbose = FALSE) {
  query <- 0
  gap <-
    gapCreation(Liste, rate = taille_gap / length(Liste), begin = start_gap)$output_vector
  cond <- 0
  
  index <- 1:length(Liste)
  
  df <- data.frame(index = index)
  
  df$Ref <- Liste
  
  score <- NULL
  
  for (i in 6:283) {
    df$signal <- compute.DTWBI_QF(
      sig = gap,
      acceptedHole = 30,
      smallHole = 3,
      S_Query = i,
      verbose = FALSE
    )
    print(length(Liste))
    print(length(df$signal))
    valeur <- compute.sim(Liste, df$signal)
    if (verbose == TRUE && i %% 20 == 0) {
      print(i)
    }
    score <- c(score, valeur)
    if (cond < valeur) {
      query = i
      cond = valeur
      if (verbose == TRUE) {
        print(i)
        print("TOP")
        print(valeur)
        
      }
    }
  }
  res <- list(var1 = query, var2 = score)
  return(res)
}
#-------------------------------------------------------------------------------

bob <- f2(tsNH4Complete, 30/length(tsNH4Complete), 3630, verbose = TRUE)

NH4 <- tsNH4Complete
si <- gapCreation(NH4, 30 / length(NH4), 3630)
for (i in 6:183) {
  print(i)
  res <- compute.DTWBI_QF(si$output_vector,30,3,i,verbose = FALSE)
  print(compute.indicateurComp(NH4,res))
}


a <- 6
b <- 183
df <- data.frame(Query = a:b,
                 R2 = NA,
                 FA2 = NA,
                 FB = NA,
                 FS = NA,
                 NMSE = NA,
                 NMSE_O = NA,
                 MG = NA,
                 VG = NA,
                 simAire = NA,
                 distMax = NA,
                 all = NA)

for (i in a:b) {
  print(i)
  DTWBI <- compute.DTWBI_QF(
    si$output_vector,
    acceptedHole = 30,
    smallHole = 0,
    S_Query = i,
    verbose = FALSE
  )
  temp <- compute.indicateurComp(NH4, DTWBI, v = F)
  # print(temp)
  for (j in 1:length(temp)+1) {
    df[i-(a-1),j] <- temp[[j-1]]
  }
}
#-------------------------------------------------------------------------------
NH4_Complete <- tsNH4Complete

gap <- gapCreation(tsNH4Complete, rate = 30/4552, begin = 3630)$output_vector

index <- 1:length(NH4_Complete)

df <- data.frame(index = index)

df$NH4_Complete <- NH4_Complete

df$signal2 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                               smallHole = 3, S_Query = 2, verbose = FALSE)

df$signal51 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                smallHole = 3, S_Query = 51, verbose = FALSE)

df$signal60 <- compute.DTWBI_QF(sig = gap, acceptedHole = 30,
                                smallHole = 3, S_Query = 60, verbose = FALSE)

ggplot(df, aes(x = index)) +
  geom_line(aes(y = NH4_Complete, color = "Complete"), size = 2) +
  geom_line(aes(y = signal2, color = "Query = 2")) +
  geom_line(aes(y = signal51, color = "Query = 51")) +
  geom_line(aes(y = signal60, color = "Query = 60")) +
  scale_x_continuous(limits = c(index[3625], index[3665])) +
  scale_y_continuous(limits = c(5, 17.5), expand = c(0, 0))