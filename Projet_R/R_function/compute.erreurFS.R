#---------------------------------------------------------------------
#' Calcule FS fraction de la variance. fractional variance
#' permet de savoir si un modèle est acceptable. Pour cela, le FS doit se rapprocher de zéro. 
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FS
compute.erreurFS<-function(courbe1, courbe2,verbose=F){
  sd1=sd(courbe1);
  sd2=sd(courbe2);
  var1=sd1^2;
  var2=sd2^2;
  FS <- 2*(var1-var2)/(var1+var2) 
  if (verbose){
    if(abs(FS)<0.5) {print("acceptable model");
    }else{print("non acceptable FS");}
  }
  out<- FS
}