#---------------------------------------------------------------------
#' MG : Calcule biais de l'erreur geometrique Ahuja and Kumar (1996) 
#' exp(ln meanc1-ln mean c2)
#' acceptable si 0.75<MG< 1.25
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2 predite  de même taille
#' @return erreur
compute.erreurMG<-function(courbe1, courbe2, verbose=F){
  
  N <- length(courbe1);
  courbe1[courbe1==0]=0.000001;
  courbe2[courbe2==0]=0.000001;
  m1 <- mean(log(courbe1));
  m2 <- mean(log(courbe2));
  res <- exp(m1-m2);
  if (verbose){
    if(res<=1.25 & res>=0.75) {print("acceptable MG error");
    }else{print("non acceptable MG error");}
  }
  out<- res
}