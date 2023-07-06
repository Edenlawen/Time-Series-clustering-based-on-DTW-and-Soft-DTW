#---------------------------------------------------------------------
#' Calcule similarite entre deux courbes de même taille basees sur l'inverse de l'aire entre les deux courbes 
#' sim= 1/N sum 1/1+(x-y)^2
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 
#' @param Courbe2  de même taille
#' @return similarite
compute.similariteAire<-function(courbe1, courbe2){
  
  simi <- 0;
  simi <- sum(1/(1+(courbe1-courbe2)^2));
  simi <- simi/length(courbe1)
}


