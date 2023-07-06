#---------------------------------------------------------------------
#' Calcule FA2 pourcentage de données 
#' représente la fraction (pourcentage) de données qui satisfait 0,5 ≤ (C2/C1) <= 2 
#' Un modèle est considéré comme parfait lorsque son FA2 est proche de 1.
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FA2
compute.erreurFA2<-function(courbe1, courbe2,verbose=F){
  ratio=(courbe2/courbe1);
  fraction=ratio[ratio>=0.5 & ratio<=2];
  FA2 <- length(fraction)/length(courbe1) 
  if (verbose){
    if(FA2>0.8) {print("good model");
    }else{print("important number of different point");}
  }
  out<- FA2
}