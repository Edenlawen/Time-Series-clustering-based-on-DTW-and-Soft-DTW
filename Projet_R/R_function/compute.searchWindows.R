#' This program allows you to search for data windows similar to an input or specified data set.
#' Author Merlin Roudier v17/07/2023
#' @param
#' @return

compute.searchWindowsWithGap <-
  function(data,
           gapTaille,
           gapStart,
           queryTaille,
           verbose = F) {
    print("Chercher fenetre viable")
    # gapTaille <- 7
    # gapStart <- 400
    dataModif <-
      gapCreation(data, gapTaille / length(data), gapStart)$output_vector
    # queryTaille <- 12
    featureRef <-
      globalfeatures(dataModif[(gapStart - queryTaille):(gapStart - 1)])
    queryRef <- dataModif[(gapStart - queryTaille):(gapStart - 1)]
    fenetresViable <- data.frame("queryRef" = queryRef)
    repRef <- data[gapStart:(gapStart + gapTaille - 1)]
    reponseViable <- data.frame("repRef" = repRef)
    debut <- 1
    fin <- debut + queryTaille
    
    if (length(data) < 1000) {
      step_threshold = 2
      stepBase <- 2
    } else{
      if (length(data) > 10000) {
        step_threshold = 50
        stepBase <- 50
      } else{
        step_threshold = 10
        stepBase <- 10
      }
    }
    
    if (length(data) < 10000) {
      threshold_cos = 0.9995
    } else {
      threshold_cos = 0.995
    }
    
    while ((fin + queryTaille + gapTaille) < length(data)) {
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
          # print(debut)
          # print(cosCompare)
          fenetresViable <- cbind(fenetresViable,
                                  queryTemp)
          colnames(fenetresViable)[ncol(fenetresViable)] <-
            paste0("Debut = ", debut)
          
          repTemp <- dataModif[fin:(fin + gapTaille - 1)]
          reponseViable <- cbind(reponseViable, repTemp)
          colnames(reponseViable)[ncol(reponseViable)] <-
            paste0("Debut = ", debut)
        }
      }
      debut <- debut + step_threshold
      fin <- fin + step_threshold
      
      if (fin >= length(data)) {
        if (length(fenetresViable) < 200) {
          if (step_threshold == 1) {
            threshold_cos <- 0.95
            step_threshold <- stepBase
          }
          debut <- 1
          fin <- debut + queryTaille
          step_threshold <-
            step_threshold - floor(step_threshold / 2)
          print(
            paste0(
              "Pas assez de fenetre viable, votre step_threshold est passé à ",
              step_threshold
            )
          )
          rm(fenetresViable, reponseViable)
          fenetresViable <- data.frame("queryRef" = queryRef)
          reponseViable <- data.frame("repRef" = repRef)
        }
      }
    }
    fenetresViable <- subset(fenetresViable, select = -1)
    reponseViable <- subset(reponseViable, select = -1)
    
  }