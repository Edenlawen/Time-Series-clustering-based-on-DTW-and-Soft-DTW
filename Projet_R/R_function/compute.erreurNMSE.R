#---------------------------------------------------------------------
#' Calcule erreur quadratique normalisee NMSE 
#' Elle met en évidence la dispersion (l’éparpillement) de données. (c1-c2)^2/meanC1.meanC2
#' acceptable si erreur<0,4
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 
#' @param Courbe2  de même taille
#' @return erreur
compute.erreurNMSE<-function(courbe1,courbe2 ,verbose=F){
  
  N <- length(courbe1);
  Dif <- 0; m1 <- 0; m2 <- 0;
  Dif <- mean((courbe1-courbe2)^2);
  m1 <- mean(courbe1);
  m2 <- mean(courbe2);
  res <- Dif/(m1*m2)
  if (verbose){
    if(res<0.4) {print("acceptable NMSE error");
    }else{print("non acceptable NMSE error");}
  }
  out<- res
}