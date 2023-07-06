#---------------------------------------------------------------------
#' VG : Calcule variance de l'erreur geometrique Ahuja and Kumar (1996) 
#' VG=exp(mean(ln c1-ln c2)^2))
#' acceptable si geometric mean variance
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2 predite  de même taille
#' @return erreur
compute.erreurVG<-function(courbe1, courbe2,verbose=F){
  
  N <- length(courbe1);
  courbe1[courbe1==0]=0.000001;
  courbe2[courbe2==0]=0.000001;
  m1 <- log(courbe1);
  m2 <- log(courbe2);
  res <- exp(mean((m1-m2)^2));
  if (verbose){
    if(res<=1.25 & res>=0.75) {print("acceptable VG error");
    }else{print("non acceptable VG error");}
  }
  out<- res
}