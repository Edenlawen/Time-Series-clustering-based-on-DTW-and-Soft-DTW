#---------------------------------------------------------------------
#' Calcule biais fractionnel FB 
#' Ce paramètre détermine si les valeurs prédites sont surestimées ou sous-estimées par rapport à celles observées. 
# Un modèle est jugé parfait quand son FB tend vers zéro. acceptable lorsque −0,3 ≤ FB ≤ 0,3
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FB
compute.erreurFB<-function(courbe1, courbe2,verbose=F){
  m1=mean(courbe1);
  m2=mean(courbe2);
  FB <- 2*(m1-m2)/(m1+m2) 
  if (verbose){
    if(abs(FB)<0.3) {print("acceptable FB");
    }else{print("non acceptable FB");}
  }
  out<- FB
}