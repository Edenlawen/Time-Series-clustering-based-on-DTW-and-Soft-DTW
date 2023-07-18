#ensemble d'indicateurs de comparaison de courbes observees/predites


#---------------------------------------------------------------------
#' Calcule similarite entre deux courbes de même taille basees sur l'inverse
#' de l'aire entre les deux courbes
#' sim= 1/N sum 1/1+(x-y)^2
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1
#' @param Courbe2  de même taille
#' @return similarite
compute.similariteAire <- function(courbe1, courbe2) {
  simi <- 0
  
  simi <- sum(1 / (1 + (courbe1 - courbe2) ^ 2))
  
  simi <- simi / length(courbe1)
}


#---------------------------------------------------------------------
#' Calcule R coefficient de correlation entre les courbes
#' mesure la ressemblance entre les données observées et les données prédites.
#' Un modèle est considéré comme acceptable lorsque R est proche de 1.
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return r2
compute.Rcorrelation <- function(courbe1, courbe2, verbose = F) {
  R = cor.test(courbe2, courbe1)
  
  if (verbose) {
    if (R$estimate ^ 2 > 0.9 &
        R$p.value < 0.05) {
      print("acceptable model")
      
    } else{
      print("not acceptable model R^2<0.9")
    }
  }
  res = R$estimate ^ 2
  out <- res
}


#---------------------------------------------------------------------
#' Calcule distance maximale entre les deux courbes
#' @Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1
#' @param Courbe2  de même taille
#' @return dmax
compute.distMaxi <- function(courbe1, courbe2) {
  dMax <- 0
  
  dMax <- max(abs(courbe2 - courbe1))
}

#---------------------------------------------------------------------
#' Calcule FA2 pourcentage de données
#' représente la fraction (pourcentage) de données qui satisfait 0,5 ≤ (C2/C1) <= 2
#' Un modèle est considéré comme parfait lorsque son FA2 est proche de 1.
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FA2
compute.erreurFA2 <- function(courbe1, courbe2, verbose = F) {
  ratio = (courbe2 / courbe1)
  
  fraction = ratio[ratio >= 0.5 & ratio <= 2]
  
  FA2 <- length(fraction) / length(courbe1)
  if (verbose) {
    if (FA2 > 0.8) {
      print("good model")
      
    } else{
      print("important number of different point")
    }
  }
  out <- FA2
}

#---------------------------------------------------------------------
#' Calcule biais fractionnel FB
#' Ce paramètre détermine si les valeurs prédites sont surestimées
#' ou sous-estimées par rapport à celles observées.
#' Un modèle est jugé parfait quand son FB tend vers zéro. acceptable
#' lorsque −0,3 ≤ FB ≤ 0,3
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FB
compute.erreurFB <- function(courbe1, courbe2, verbose = F) {
  m1 = mean(courbe1)
  
  m2 = mean(courbe2)
  
  FB <- 2 * (m1 - m2) / (m1 + m2)
  if (verbose) {
    if (abs(FB) < 0.3) {
      print("acceptable FB")
      
    } else{
      print("non acceptable FB")
    }
  }
  out <- FB
}

#---------------------------------------------------------------------
#' Calcule FS fraction de la variance. fractional variance
#' permet de savoir si un modèle est acceptable.
#' Pour cela, le FS doit se rapprocher de zéro.
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2  predite
#' @return FS
compute.erreurFS <- function(courbe1, courbe2, verbose = F) {
  sd1 = sd(courbe1)
  
  sd2 = sd(courbe2)
  
  var1 = sd1 ^ 2
  
  var2 = sd2 ^ 2
  
  FS <- 2 * (var1 - var2) / (var1 + var2)
  if (verbose) {
    if (abs(FS) < 0.5) {
      print("acceptable model")
      
    } else{
      print("non acceptable FS")
    }
  }
  out <- FS
}

#---------------------------------------------------------------------
#' MG : Calcule biais de l'erreur geometrique Ahuja and Kumar (1996)
#' exp(ln meanc1-ln mean c2)
#' acceptable si 0.75<MG< 1.25
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2 predite  de même taille
#' @return erreur
compute.erreurMG <- function(courbe1, courbe2, verbose = F) {
  N <- length(courbe1)
  
  courbe1[courbe1 == 0] = 0.000001
  
  courbe2[courbe2 == 0] = 0.000001
  
  m1 <- mean(log(courbe1))
  
  m2 <- mean(log(courbe2))
  
  res <- exp(m1 - m2)
  
  if (verbose) {
    if (res <= 1.25 & res >= 0.75) {
      print("acceptable MG error")
      
    } else{
      print("non acceptable MG error")
    }
  }
  out <- res
}

#---------------------------------------------------------------------
#' Calcule erreur quadratique normalisee NMSE
#' Elle met en évidence la dispersion (l’éparpillement) de données.
#' (c1-c2)^2/meanC1.meanC2
#' acceptable si erreur<0,4
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1
#' @param Courbe2  de même taille
#' @return erreur
compute.erreurNMSE <- function(courbe1, courbe2 , verbose = F) {
  N <- length(courbe1)
  
  Dif <- 0
  m1 <- 0
  m2 <- 0
  
  Dif <- mean((courbe1 - courbe2) ^ 2)
  
  m1 <- mean(courbe1)
  
  m2 <- mean(courbe2)
  
  res <- Dif / (m1 * m2)
  if (verbose) {
    if (res < 0.4) {
      print("acceptable NMSE error")
      
    } else{
      print("non acceptable NMSE error")
    }
  }
  out <- res
}

#---------------------------------------------------------------------
#' Calcule erreur quadratique normalisee ratio des aires erreur/original
#' sim= sqrt(1/N sum (x-y)^2 / 1/N sum (x)^2)
#' satisfaisant si tend vers 0 mauvaise lorsque tend vers 1 tres mauvaise
#' superieur à  acceptable si erreur<0,4
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1
#' @param Courbe2  de même taille
#' @return erreur
compute.erreurqn <- function(courbe1, courbe2, verbose = F) {
  N <- length(courbe1)
  
  aire1 <- 0
  
  aire1 <- sum((courbe1 - courbe2) ^ 2)
  
  aire2 <- 0
  
  aire2 <- sum((courbe1) ^ 2)
  
  simi <- sqrt(aire1 / aire2)
  
  if (verbose) {
    if (simi < 0.4) {
      print("acceptable NMSE_O error")
      
    } else{
      print("non acceptable NMSE_O error")
    }
  }
  out <- simi
  
  
}


#---------------------------------------------------------------------
#' VG : Calcule variance de l'erreur geometrique Ahuja and Kumar (1996)
#' VG=exp(mean(ln c1-ln c2)^2))
#' acceptable si geometric mean variance
#' Author Emilie Poisson Caillault v22/05/2015
#' @param courbe1 observee
#' @param Courbe2 predite  de même taille
#' @return erreur
compute.erreurVG <- function(courbe1, courbe2, verbose = F) {
  N <- length(courbe1)
  
  courbe1[courbe1 == 0] = 0.000001
  
  courbe2[courbe2 == 0] = 0.000001
  
  m1 <- log(courbe1)
  
  m2 <- log(courbe2)
  
  res <- exp(mean((m1 - m2) ^ 2))
  
  if (verbose) {
    if (res <= 1.25 & res >= 0.75) {
      print("acceptable VG error")
      
    } else{
      print("non acceptable VG error")
    }
  }
  out <- res
}




#---------------------------------------------------------------------
#' Calcul de différents indicateurs de comparaison pour des courbes de même taille
#' Author Emilie Poisson Caillault v22/05/2015 | Modifier par Merlin
#' @param courbe1 observée
#' @param Courbe2  prédite
#' @return indicateurs
compute.indicateurComp <-
  function(courbe1,
           courbe2,
           par.R2 = 0.9,
           par.FA2 = 0.8,
           par.FB = 0.3,
           par.FS = 0.05,
           par.NMSE = 0.4,
           par.simAire = 0.95,
           par.MGmin = 0.75,
           par.MGmax = 1.25,
           v = F) {
    mini <- min(courbe1, courbe2)
    if (mini < 0) {
      courbe1 <- courbe1 + abs(mini)
      courbe2 <- courbe2 + abs(mini)
    }
    res = NULL
    
    res = c("R2" = compute.Rcorrelation(courbe1, courbe2, verbose = v))
    
    res = c(res, "FA2" = compute.erreurFA2(courbe1, courbe2, verbose = v))
    
    res = c(res, "FB" = compute.erreurFB(courbe1, courbe2, verbose = v))
    
    res = c(res, "FS" = compute.erreurFS(courbe1, courbe2, verbose = v))
    
    res = c(res, "NMSE" = compute.erreurNMSE(courbe1, courbe2, verbose = v))
    
    res = c(res, "NMSE_O" = compute.erreurqn(courbe1, courbe2, verbose = v))
    
    res = c(res, "MG" = compute.erreurMG(courbe1, courbe2, verbose = v))
    
    res = c(res, "VG" = compute.erreurVG(courbe1, courbe2, verbose = v))
    
    res = c(res, "simAire" = compute.similariteAire(courbe1, courbe2))
    
    res = c(res, "distMax" = compute.distMaxi(courbe1, courbe2))
    
    test = T
    
    test = test &
      res["R2.cor"] >= par.R2 &
      res["FA2"] >= par.FA2 &
      abs(res["FB"]) < par.FB &
      abs(res["FS"]) < par.FS & res["NMSE"] < par.NMSE
    
    test = test &
      res["simAire"] > par.simAire &
      res["MG"] >= par.MGmin &
      res["VG"] >= par.MGmin  &
      res["VG"] <= par.MGmax
    
    totT <- 0
    totF <- 0
    if (res["R2.cor"] >= par.R2) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["FA2"] >= par.FA2) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (abs(res["FB"]) < par.FB) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (abs(res["FS"]) < par.FS) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["NMSE"] < par.NMSE) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["simAire"] > par.simAire) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["MG"] >= par.MGmin) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["VG"] >= par.MGmin) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    if (res["VG"] <= par.MGmax) {
      totT <- totT + 1
    } else{
      totF <- totF + 1
    }
    
    if (test) {
      res = c(res, "all" = 1)
      
    } else{
      res = c(res, "all" = -1)
      
    }
    
    res = c(res, "Nombre de conditions vraies" = totT)
    res = c(res, "Nombre de conditions fausses" = totF)
    
    out <- res
    
  }