#---------------------------------------------------------------------
#' Calcule R coefficient de correlation entre les courbes 
#'  mesure la ressemblance entre les données observées et les données prédites.
#' Un modèle est considéré comme acceptable lorsque R est proche de 1.  
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return r2
compute.Rcorrelation<-function(courbe1, courbe2, verbose=F){
  R=cor.test(courbe2,courbe1);
  if (verbose){
    if(R$estimate^2>0.9 & R$p.value<0.05) {print("acceptable model");
    }else{print("not acceptable model R^2<0.9");}
  }
  res=R$estimate^2
  out<- res
}