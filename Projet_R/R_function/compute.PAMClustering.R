library(cluster)
library(factoextra)

compute.PAMClustering <-
  function(data,
           matrice,
           graphic = TRUE,
           queryRef = NULL,
           nbcluster = NULL) {
    if (is.null(nbcluster)) {
      nbcluster <-
        fviz_nbclust(t(data),
                     pam,
                     diss = matrice,
                     method = "silhouette",)
      nbcluster <- nbcluster$data$y
      nbcluster <- which.max(nbcluster)
    }
    
    resultatPam <-
      pam(matrice,
          nbcluster,
          diss = TRUE,
          cluster.only = TRUE)
    
    if (graphic) {
      queryTaille <- length(fenetresViable[, 1])
      # Partie plot de courbe comme Q1; mediane; Q3
      dfPAM <- list()
      for (i in 1:nbcluster) {
        dfPAM <-
          append(dfPAM, list(data.frame("index" = 1:queryTaille)), after = length(dfPAM))
      }
      for (i in 1:length(data)) {
        df_selected <- dfPAM[[resultatPam[i]]]
        df_selected <- cbind(df_selected, data[i])
        dfPAM[[resultatPam[i]]] <- df_selected
      }
      
      for (i in 1:nbcluster) {
        df_selected <- dfPAM[[i]]
        df_selected <- t(df_selected)
        df_selected <- df_selected[2:nrow(df_selected),]
        dfPAM[[i]] <- df_selected
      }
      
      summaryDTW <- list()
      for (i in 1:nbcluster) {
        summaryDTW <-
          append(summaryDTW, list(data.frame("index" = 1:queryTaille)), after = length(summaryDTW))
      }
      
      
      for (i in 1:nbcluster) {
        Q1 <- NULL
        med <- NULL
        Q3 <- NULL
        df_selected <- summaryDTW[[i]]
        df_clust <- dfPAM[[i]]
        if (length(df_clust) != queryTaille) {
          for (j in 1:queryTaille) {
            Q1 <- c(Q1, quantile(df_clust[, j], 0.25))
            med <- c(med, quantile(df_clust[, j], 0.5))
            Q3 <- c(Q3, quantile(df_clust[, j], 0.75))
          }
          df_selected <- cbind(df_selected, queryRef)
          df_selected <- cbind(df_selected, Q1)
          df_selected <- cbind(df_selected, med)
          df_selected <- cbind(df_selected, Q3)
          summaryDTW[[i]] <- df_selected
        } else{
          Q1 <- df_clust[1:queryTaille]
          med <- df_clust[1:queryTaille]
          Q3 <- df_clust[1:queryTaille]
          df_selected <-
            cbind(df_selected, queryRef)
          df_selected <-
            cbind(df_selected, Q1)
          df_selected <-
            cbind(df_selected, med)
          df_selected <-
            cbind(df_selected, Q3)
          summaryDTW[[i]] <- df_selected
        }
      }
      
      
      plotPAMDTW <- list()
      for (i in 1:nbcluster) {
        t <- summaryDTW[[i]]
        t <- gather(t, key = variable, value = value,-index)
        p <-
          ggplot(data = t, aes(x = index, y = value, color = variable))
        p <-
          p + geom_line() + labs(x = "Jour",
                                 y = "Valeur",
                                 title = paste("PAM, cluster"))
        
        # Affichage du graphique avec plotly
        plotPAMDTW[[i]] <- ggplotly(p)
      }
      
      print(subplot(plotPAMDTW, nrows = 2))
    }
    
    return(list(nbcluster, resultatPam))
  }