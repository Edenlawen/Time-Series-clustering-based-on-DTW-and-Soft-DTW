#' @title Estimating global features of a univariate signal
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault v2019.04.18
#' @description Computes global features of a univariate signal, used as input for threshold and window definition in DTWBI algorithm.
#'  Features computed are:
#'  \itemize{
#'  \item{minx: }{minimum value of the input vector}
#'  \item{maxx: }{maximum value of the input vector}
#'  \item{avg: }{average value of the input vector}
#'  \item{medianx: }{median of the input vector}
#'  \item{sdx: }{standard deviation of the input vector}
#'  \item{momx: }{skewness}
#'  \item{nop: }{number of peaks of the input vector}
#'  \item{len: }{length of the input vector}
#'  \item{entro: }{measure of entropy}
#'  }
#' @param X signal
#' @return returns a matrix with one row, each column giving the value of corresponding estimated feature.
#' @importFrom entropy entropy
#' @importFrom e1071 skewness
#' @keywords internal
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1]
#' gf <- .globalfeatures(X)

.findPeaks <-function(x, thresh = 0) {
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  if (!missing(thresh)) {
    (pks[x[pks] - x[pks + 1] > thresh]) || (pks[x[pks] - x[pks - 1] > thresh])
  }
  else pks
}

globalfeatures<-function(X){
  minx <- min(X,na.rm=T); 
  maxx <- max(X,na.rm=T); 
  avg <- mean(X,na.rm=T) ;
  medianx <- median(X,na.rm=T); 
  sdx <- sd(X,na.rm=T); 
  mom3 <- e1071::skewness(X,na.rm = T); 
  nop <- length(.findPeaks(X)); 
  len <- length(X)
  entro <- entropy::entropy(as.vector(table(X)), method="ML"); 
  out <- c(minx,maxx,avg,medianx,sdx,mom3,nop,len,entro)
  # out <- format(out, digits=2, nsmall=2)
  return(out)
}