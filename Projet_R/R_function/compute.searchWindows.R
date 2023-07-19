#' This program allows you to search for data windows similar to an input or specified data set.
#' Author Merlin Roudier v.../07/2023
#' @param
#' @return

compute.searchWindows <-
  function(data,
           query,
           verbose = F) {
    queryTaille <- len(query)
    featureRef <-
      globalfeatures(query)
    fenetresViable <- data.frame("query" = query)
    debut <- 1
    fin <- debut + queryTaille
    
    if (length(donnee) < 1000) {
      borne_a <- 1
      borne_b <- 2
      step_threshold <- sample(borne_a:borne_b, 1)
    } else{
      if (length(donnee) > 10000) {
        borne_a <- 11
        borne_b <- 50
        step_threshold <- sample(borne_a:borne_b, 1)
      } else{
        borne_a <- 3
        borne_b <- 10
        step_threshold <- sample(borne_a:borne_b, 1)
      }
    }
    
    threshold_cos <- 0.95
    puiss <- 3
    
    while (((fin + queryTaille + gapTaille) < length(donnee))) {
      if (!(
        debut %in% seq(
          gapStart - queryTaille - queryTaille,
          gapStart + gapTaille + queryTaille
        )
      )) {
        featureTemp <- globalfeatures(dataModif[debut:(fin - 1)])
        queryTemp <- dataModif[debut:(fin - 1)]
        cosCompare <- abs(cosine(featureRef, featureTemp))
        
        if (!is.na(cosCompare) && cosCompare >= threshold_cos) {
          fenetresViable <- cbind(fenetresViable,
                                  queryTemp)
          colnames(fenetresViable)[ncol(fenetresViable)] <-
            paste0("Debut = ", debut)
        }
      }
      step_threshold <- sample(borne_a:borne_b, 1)
      debut <- debut + step_threshold
      fin <- fin + step_threshold
      
      if ((fin + queryTaille + gapTaille) >= length(donnee)) {
        if (length(fenetresViable) > 150) {
          debut <- 1
          fin <- debut + queryTaille
          threshold_cos <- threshold_cos + (45 * (10 ^ (-puiss)))
          puiss <- puiss + 1
          if (verbose) {
            print(threshold_cos)
          }
          rm(fenetresViable, reponseViable)
          fenetresViable <- data.frame("query" = query)
        }
      }
    }
    fenetresViable <- subset(fenetresViable, select = -1)
  }