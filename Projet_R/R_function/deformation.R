library(splines)
library(stats)
#-------------------------------------------------------------------------------
deformation_polynomiale <- function(courbe_initiale, courbe_cible, degre_polynome) {
  # Ajustement d'un polynôme de degré donné
  modele_polynome <- lm(courbe_cible ~ poly(courbe_initiale, degre_polynome, raw = TRUE))
  
  # Prédiction de la courbe déformée
  courbe_deformee <- predict(modele_polynome, newdata = data.frame(courbe_initiale))
  
  return(courbe_deformee)
}



#-------------------------------------------------------------------------------
deformation_par_deplacement <-
  function(courbe_initiale, courbe_cible) {
    nb_points <- length(courbe_initiale)
    courbe_deformee <- numeric(nb_points)
    
    for (i in 1:nb_points) {
      courbe_deformee[i] <-
        optim(courbe_initiale[i], function(x)
          (x - courbe_cible[i]) ^ 2)$par
    }
    
    return(courbe_deformee)
  }


#-------------------------------------------------------------------------------
deformation_par_warping <- function(courbe_initiale, courbe_cible) {
  nb_points <- length(courbe_initiale)
  courbe_deformee <- numeric(nb_points)
  
  for (i in 1:nb_points) {
    courbe_deformee[i] <-
      courbe_initiale[i] + (courbe_cible[i] - courbe_initiale[i])
  }
  
  return(courbe_deformee)
}


#-------------------------------------------------------------------------------
deformation_par_optimisation <-
  function(courbe_initiale, courbe_cible) {
    fonction_objectif <-
      function(parametres)
        sum((parametres * courbe_initiale - courbe_cible) ^ 2)
    resultats_optim <- optim(par = 1, fn = fonction_objectif)$par
    courbe_deformee <- resultats_optim * courbe_initiale
    
    return(courbe_deformee)
  }
