#' This program allows you to search for data windows similar to an input or specified data set.
#' Author Merlin Roudier v21/07/2023
#' @param data
#' @param queryRef
#' @param nbWindows
#' @param seeds
#' @param random
#' @param verbose
#' @return
source("R_function/globalF.r")
source("R_function/EC_compareCourbe.r")
library(lsa)

# compute.searchWindows <-
#   function(data,
#            queryRef,
#            startqueryRef,
#            nbWindows = 500,
#            seeds = NULL,
#            random = TRUE,
#            verbose = F) {
#     queryRefTaille <- length(queryRef)
#     fenetresViable <- data.frame("queryRefRef" = queryRef)
#     featureRef <- globalfeatures(queryRef)
#     debut <- 1
#     fin <- debut + queryRefTaille
#     
#     if (random) {
#       set.seed(seeds)
#       if (length(data) < 1000) {
#         borne_a <- 1
#         borne_b <- 2
#         step_threshold <- sample(borne_a:borne_b, 1)
#       } else{
#         if (length(data) > 10000) {
#           borne_a <- 11
#           borne_b <- 50
#           step_threshold <- sample(borne_a:borne_b, 1)
#         } else{
#           borne_a <- 3
#           borne_b <- 10
#           step_threshold <- sample(borne_a:borne_b, 1)
#         }
#       }
#     } else{
#       if (length(data) < 1000) {
#         step_threshold <- 1
#       } else{
#         if (length(data) > 10000) {
#           step_threshold <- 50
#         } else{
#           step_threshold <- 10
#         }
#       }
#     }
#     
#     threshold_cos <- 0.95
#     puiss <- 3
#     step_threshold <- 2
#     
#     cos_score <- c()
#     deb_vect <- c()
#     
#     while (((fin + queryRefTaille) < length(data))) {
#       if (!(debut %in% seq(startqueryRef,
#                            startqueryRef + queryRefTaille))) {
#         featureTemp <- globalfeatures(data[debut:(fin - 1)])
#         cosCompare <- abs(cosine(featureRef, featureTemp))
#         
#         if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
#           cos_score <- c(cos_score, cosCompare)
#           deb_vect <- c(deb_vect, debut)
#         }
#       }
#       
#       if (random) {
#         step_threshold <- sample(borne_a:borne_b, 1)
#       }
#       debut <- debut + step_threshold
#       fin <- debut + queryRefTaille
#       
#       if ((fin + queryRefTaille) >= length(data)) {
#         if (length(cos_score) > nbWindows &&
#             length(cos_score) < (nbWindows + 2000)) {
#           while (length(cos_score) > nbWindows &&
#                  length(cos_score) < (nbWindows + 2000)) {
#             threshold_cos <- threshold_cos + (45 * (10 ^ (-puiss)))
#             puiss <- puiss + 1
#             if (verbose) {
#               print(threshold_cos)
#             }
#             critere <- function(x, y) {
#               x < y
#             }
#             cos_score <-
#               cos_score[!critere(cos_score, threshold_cos)]
#             deb_vect <- deb_vect[!critere(cos_score, threshold_cos)]
#           }
#         }
#       }
#     }
#     for (i in 1:length(deb_vect)) {
#       debut <- deb_vect[i]
#       fin <- debut + queryRefTaille
#       if (compute.indicateurComp(queryref, data[debut:(fin - 1)])[12] >= 7) {
#         fenetresViable <- cbind(fenetresViable,
#                                 queryRefTemp <-
#                                   data[debut:(fin - 1)])
#         colnames(fenetresViable)[ncol(fenetresViable)] <-
#           paste0("Debut = ", deb_vect[i])
#       }
#     }
#     fenetresViable <- subset(fenetresViable, select = -1)
#     return(fenetresViable)
#   }


compute.searchWindows <-
  function(data,
           queryRef,
           startqueryRef,
           nbWindows = 500,
           seeds = NULL,
           random = TRUE,
           verbose = F) {
    queryRefTaille <- length(queryRef)
    fenetresViable <- data.frame("queryRefRef" = queryRef)
    featureRef <- globalfeatures(queryRef)
    debut <- 1
    fin <- debut + queryRefTaille
    
    if (random) {
      set.seed(seeds)
      if (length(data) < 1000) {
        borne_a <- 1
        borne_b <- 2
        step_threshold <- sample(borne_a:borne_b, 1)
      } else{
        if (length(data) > 10000) {
          borne_a <- 11
          borne_b <- 50
          step_threshold <- sample(borne_a:borne_b, 1)
        } else{
          borne_a <- 3
          borne_b <- 10
          step_threshold <- sample(borne_a:borne_b, 1)
        }
      }
    } else{
      if (length(data) < 10000) {
        step_threshold <- 1
      } else{
        if (length(data) > 1000000) {
          step_threshold <- 10
        } else{
          step_threshold <- 5
        }
      }
    }
    
    threshold_cos <- 0.95
    puiss <- 3
    
    cos_score <- c()
    deb_vect <- c()
    
    while (((fin + queryRefTaille) < length(data))) {
      if (!(debut %in% seq(startqueryRef,
                           startqueryRef + queryRefTaille))) {
        featureTemp <- globalfeatures(data[debut:(fin - 1)])
        cosCompare <- abs(cosine(featureRef, featureTemp))
        
        if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
          cos_score <- c(cos_score, cosCompare)
          deb_vect <- c(deb_vect, debut)
        }
      }
      
      if (random) {
        step_threshold <- sample(borne_a:borne_b, 1)
      }
      debut <- debut + step_threshold
      fin <- debut + queryRefTaille
    }
    for (i in 1:length(deb_vect)) {
      debut <- deb_vect[i]
      fin <- debut + queryRefTaille
      if (compute.indicateurComp(queryref, data[debut:(fin - 1)])[12] >= 7) {
        fenetresViable <- cbind(fenetresViable,
                                queryRefTemp <-
                                  data[debut:(fin - 1)])
        colnames(fenetresViable)[ncol(fenetresViable)] <-
          paste0("Debut = ", deb_vect[i])
      }
    }
    fenetresViable <- subset(fenetresViable, select = -1)
    return(fenetresViable)
  }