#---------------------------------------------------------------------
#' Calcule erreur quadratique normalisee ratio des aires erreur/original 
#' sim= sqrt(1/N sum (x-y)^2 / 1/N sum (x)^2)
#' satisfaisant si tend vers 0 mauvaise lorsque tend vers 1 tres mauvaise superieur à 
#' acceptable si erreur<0,4
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 
#' @param Courbe2  de même taille
#' @return erreur
compute.erreurqn<-function(courbe1, courbe2,verbose=F){
  
  N <- length(courbe1);
  aire1 <- 0;
  aire1 <- sum((courbe1-courbe2)^2);
  aire2 <- 0;
  aire2 <- sum((courbe1)^2);
  simi <- sqrt(aire1/aire2);
  if (verbose){
    if(simi<0.4) {print("acceptable NMSE_O error");
    }else{print("non acceptable NMSE_O error");}
  }
  out<- simi
}
