#---------------------------------------------------------------------
#' Calcule distance maximale entre les deux courbes 
#' @Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 
#' @param Courbe2  de même taille
#' @return dmax
compute.distMaxi<-function(courbe1, courbe2){
  
  dMax <- 0;
  dMax <- max(abs(courbe2-courbe1))
}
