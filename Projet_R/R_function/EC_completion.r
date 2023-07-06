#---------------------------------------------------------------------
#' Compute size of consecutive hole for each time in time series 
#' return vector with same size. 0 if no hole other the size of hole
#' Author Emilie Poisson Caillault v25/01/2019
#' @param serie vector 
#' @return vector 
compute.sizeHole<-function(series){
  #find hole
  indice=which(is.na(series));
  N=length(indice);
  if(N>0){
    flagIndice=indice-c(-1,indice[-N]);
    sizeHole=1;
    consecutiveHole=rep(1,N);
    for(i in 1:N){
      #increment number of consecutive NA data
      if(flagIndice[i]==1){
        sizeHole=sizeHole+1;
        #update all indices concerned by same consecutive hole
        consecutiveHole[(i-(sizeHole-1)):i]=sizeHole;
      }else{
        sizeHole=1;
      }
    }
    sizeHole=rep(0,length(series));
    sizeHole[indice]=consecutiveHole;
  }else{
    sizeHole=0;
    print("no hole or no NA");
  }
  return(sizeHole);
}



#---------------------------------------------------------------------
#' indexes of missing hole use compute sizeHole 
#' return a dataframe with the position (pos) of the hole and its size  
#' Author Emilie Poisson Caillault v16/03/2015
#' @param serie vector 
#' @return dataframe if existing hole, else dataframe(pos=0,size=0)
compute.infoHole<-function(serie){
  vgap=compute.sizeHole(serie);
  index=NULL; i=0; posGap=NULL; sizeGap=NULL;
  if(sum(vgap)>0){
    indPos=which(vgap>0);
    index=data.frame(pos=indPos[1],size=vgap[indPos[1]]);
    i=1+vgap[indPos[1]];
    while(i<=length(indPos)){
      posGap=indPos[i];
      sizeGap=vgap[posGap]
      index[nrow(index)+1,]=c(posGap,sizeGap);
      i=i+sizeGap;
      posGap=NULL; sizeGap=NULL;
    }
  }else{
    index=data.frame(pos=0,size=0);
  }
  return(index);
}


#---------------------------------------------------------------------
#' hist of hole use compute sizeHole 
#' return the hist vector 
#' Author Emilie Poisson Caillault v16/03/2015
#' @param fig T or F, T=plot hist
#' @param serie vector 
#' @param verbose true or false to print some details
#' @return vector 
compute.histSizeHole<-function(series, fig=F,verbose=F){
  val=NULL;count=NULL;
  #find hole
  holeSize=compute.sizeHole(series);
  if(length(holeSize)>1){
    unique.holeSize=holeSize;
    N=length(series);
    val=NULL; count=NULL;
    for (i in 2:N){
      if(holeSize[i]>0)
        if (holeSize[i]==holeSize[i-1]) unique.holeSize[i]=-1;
    }
    #count
    maxi=max(unique.holeSize);
    mini=0;
    val=seq(mini,maxi,1);
    count=sapply(val,FUN=function(x){length(which(unique.holeSize==x))})
    if(fig){barplot(count,names.arg=val,xlab="hole size",ylab="frequency");
      if(verbose) {print(cbind(val,count))}
    }
  }else{
    if(verbose==T)
      print("no hole or no NA")
    val=0; count=length(series);
  }
  return(cbind(val,count));
}

#---------------------------------------------------------------------
#' compute information about missing data 
#' return a matrix of 4 index by parameters dataframe with the position (pos) of the hole and its size  
#' Author Emilie Poisson Caillault v25/01/2019
#' @param df vector or dataframe (parameter in column) 
#' @param v verbose true or false
#' @param fig show percent of na by column true or false
#' @return matrix 4 row index by parameters
infoNA<-function(df,v=T,fig=T){
  info=NULL;
  if(is.vector(df)) df=data.frame(df)
  info=sapply(df,FUN=function(x){
    nbNA=sum(is.na(x));
    PercentNA=0; largestGap=0;nbGap=0;infoGap=NULL;
    if(nbNA>0){
      PercentNA=100*nbNA/length(x);
      infoGap=compute.infoHole(x);
      largestGap=max(infoGap$size);
      if(infoGap$pos[1]>0) nbGap=nrow(infoGap);
    }
    return(c(nbNA,PercentNA,largestGap,nbGap))
  })
  rownames(info)=c("nbNA","PercentNA","largestGap","nbGap");
  if(v)
    if(fig) {
      plot(info['PercentNA',],main="% de NA by column",pch=3)
      abline(h=100,col="red")
    }
  if(v){print(info)}
  return(info);
}



# Toutes les fonctions de completion
# Derni?re  version avec accHole, minHole et double requette

#---------------------------------------------------------------------
#' Completion individual hole (no consecutive missing values) in time series 
#' Author Emilie Poisson Caillault v27/03/2019
#' @param serie vector 
#' @param verbose true or false to print some details 
#' @return vector
completion.singleHole<-function(series,verbose = F){
  completed=series;
  if(sum(is.na(series))==0){
    if (verbose) print("no completion");
    return(completed);
  }
  #compute size of Hole
  sizeHole=compute.sizeHole(series)
  #find individual hole
  indice=which(sizeHole==1);
  N=length(indice)
  if(N>0){
    #completion by mean
    debut=sapply(indice,function(x){max(1,x-1)})
    fin=sapply(indice,function(x){min(length(series),(x+1))})
    completed[indice]=rowMeans(cbind(series[debut],series[fin]),na.rm=TRUE);
  }
  return(completed);
}

#---------------------------------------------------------------------
#' Completion of time series by moving ponderated average according the hole size
#' Weigthed average according the place of the information, near it is higher missing data. 
#' Author Emilie Poisson Caillault v27/03/2019
#' @param serie vector
#' @param acceptedHole complete hole os size min or = to acceptedHole (default 10)
#' @param mult (default 10), if=0 standart weighting : 12..N.N..21 N=hole size else N*mult (as exponential)
#' @param confidence ponderation for NA completed value, if confidence >0 : weighted ponderation =1
#' @param maxWindow (default 10) max number of considered points to estimate hole. 
#' @param V = verbose 
#' @example t=1:500;s=sin(2*pi*t/100);s[10]=NA;s[1:2]=NA;s[480:500]=NA;s[210:240]=NA;
#' @example plot(s); sc=completion.adaptativeMovingAverage(s); points(sc,pch='+',col="red")
#' @return vector
completion.adaptativeMovingAverage<-function(series,acceptedHole=10,mult=10, confidence=1, maxWindow=10, V=FALSE){
  if(sum(is.na(series))==0) return(series);
  #fill individual hole is not already done
  completed=series;
  sizeHole=compute.infoHole(series); 
  if(sum(sizeHole[,2]==1)>0) completed=completion.singleHole(series);
  if(sum(is.na(completed))==0) return(completed);
  #compute size of Hole not filled
  sizeHole=compute.infoHole(completed);
  if(V){	
    print("hole information after single hole filling" )
    print(sizeHole)
  }
  sizeToFill=sizeHole[,2];
  positionHole=sizeHole[,1]
  indice=which(sizeToFill<=acceptedHole); # fill only missing value, no smoothing 
  N=length(indice)
  if(N==0) return(completed);
  #completion by moving average and taking account new inserted data
  for(i in indice){
    dateS.missing=positionHole[i] #start position of hole
    dateE.missing=positionHole[i]+sizeToFill[i]-1 #end position of hole
    s=sizeToFill[i];
    if(V) print(paste("----info hole i=",i, " : at ", dateS.missing, " to",dateE.missing))
    if(s==1){
      debut=max(1,dateS.missing-1)
      fin=min(length(series),dateS.missing+1)
      completed[date.missing]=mean(completed[debut:fin],na.rm=TRUE)
      hole.before=0;
    }else{
      fen=floor(min(maxWindow,length(series))/2)
      for(p in 1:s){
        dS=max(1,dateS.missing+p-1-fen);
        dE=min(dateS.missing+p-1+fen,length(series))
        ind.window=dS:dE;
        val.orig=series[dS:dE]
        val.window=completed[dS:dE];
        #create weigthed coef vector
        ponderation=rep(0,length(ind.window));
        Np=length(ponderation)
        q=(Np-fen+1):length(ponderation)
        if(mult>0){
          alpha=0.4;
          ponderation[fen:1]=alpha*(1-alpha)^(1:fen);
          ponderation=mult*ponderation/max(ponderation);
          ponderation[q]=alpha*(1-alpha)^(1:(length(q)));
          ponderation[q]=mult*ponderation[q]/max(ponderation[q]);
        }else{
          ponderation=1:length(ind.window);
          ponderation[q]=length(q):1;
        } 
        if (confidence) ponderation[is.na(val.orig)]=1;
        ponderation[is.na(val.window)]=0;
        completed[dateS.missing+p-1]=sum(ponderation*val.window,na.rm=TRUE)/sum(ponderation);
        if(V){
          print(paste("missing at ",p-1+dateS.missing))
          print("ponderation:");print(ponderation)
          print("series avant"); print(series[ind.window])
          print("series apres"); print(completed[ind.window])
        }
      }
    }
  }
  return(completed);
}

#---------------------------------------------------------------------
#' compute DTWBI algorithm on one signal
#' Author Emilie Poisson Caillault, Marco Hanocq v29/04/2022
#' @param sig 
#' @param acceptedHole size of accepted hole, to complete
#' @param smallHole size of considered small hole, completed by movering average
#' @param S_Query size of the query
#' @param thresh_cos_stop=0.8 to avoid filling
#' @param threshold_cos=0.995 for reduction computation, accepted cos between features
compute.DTWBI_QF<-function(sig, acceptedHole, smallHole, S_Query, thresh_cos_stop=0.8, threshold_cos = 0.995,verbose=T,fig=F, ...){
  sigc=NULL;
  out<-sig;
  if(verbose){
    print(paste("number of data:",length(sig)));
    print(paste("NA total:",sum(is.na(sig))));
  }
  if(sum(is.na(sig))>0){
    # filling single hole
    sigc=completion.singleHole(sig);
    # then small hole : missing consecutives values with number < smallHole
    s1=completion.adaptativeMovingAverage(sigc,acceptedHole = smallHole,mult = 10,confidence = 1,maxWindow = 10 ,V = F);
    s2=completion.adaptativeMovingAverage(sigc[length(sig):1],acceptedHole = smallHole,mult = 10,confidence = 1,maxWindow = 10,V = F);
    s2=s2[length(s2):1]
    sigc=(s1+s2)/2;
    rm(s1);rm(s2);
    #missing values after isolated hole completion
    vgap=compute.sizeHole(sigc);
    #searching position and size of each hole
    infoHole=compute.infoHole(sigc);
    infoHole=infoHole[order(infoHole[,2],decreasing=F), ]
    
    #reject Hole with size>acceptedHole
    if(verbose) print(infoHole)
    infoHole=infoHole[infoHole$size<=acceptedHole,]
    
    # acceptable Hole to fill by DTWBI methods
    sigDTWBI=sigc;
    if(nrow(infoHole)>0){
      for(i in 1:nrow(infoHole)){
        if(verbose){print( paste("% completion",floor(100*i/nrow(infoHole))))}
        posHole=infoHole$pos[i]
        sizeHole=infoHole$size[i]
        resDTWBI=NULL;
        try({
          resDTWBI=imputeDTWBI_QF(sigDTWBI, t_gap=posHole, T_gap=sizeHole, S_Query, thresh_cos_stop=thresh_cos_stop,threshold_cos = threshold_cos,verbose=F);
          sigc[posHole:(posHole+sizeHole-1)]=resDTWBI$imputation_window;
        })
      }
    }
  }
  out<-sigc
  return(out);
}



###### imputeDTWBI_QF  et fonction compl?mentaire #######################

require("entropy")
require("e1071")
require("lsa")
require("dtw")
require("forecast")

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

.globalfeatures<-function(X){
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
  out <- format(out, digits=2, nsmall=2)
  return(out)
}

#' @title Global threshold for missing data imputation
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault, Marco Hanocq
#' @date 2022/04/29
#' @description Finds a threshold for univariate missing data imputation in a univariate vector.
#' @import dtw
#' @importFrom lsa cosineplot(compute.DTWBI_QF(sp[1:length(sp)]),500,20,1)
#' @keywords internal

.DTW_threshold_global_univariate <- function(query, database, i_start, i_finish, T_gap_OG, F_half, step_threshold, threshold_cos, thresh_cos_stop, ...){
  if(sum(is.na(query))>0) return(NULL);
  # Initialization
  T_gap <- length(query)
  Cosine_threshold <- NULL;
  pos_i <- NULL;
  threshold_cos_temp <- threshold_cos
  gf_q <-as.numeric(.globalfeatures(query))
  ind_nan.q=which(is.nan(gf_q))
  while((length(Cosine_threshold)==0)&&(threshold_cos_temp>thresh_cos_stop)){
    i <- i_start
    # window browsing
    while(i<=i_finish){ 
      k <- i+T_gap-1
      ref <- database[i:k]
      if (F_half){
        j <- max(0,i-T_gap_OG)
        check_NA <- database[j:(i-1)]
      }
      else{
        j <- min(length(database),k+T_gap_OG+1)
        check_NA <- database[(k+1):j]
      }
      cos_value=0;
      if((sum(is.na(ref))==0)&&(sum(is.na(check_NA))==0)){
        #Check if the window or the answer contain NA
        gf_r <- as.numeric(.globalfeatures(ref))
        ind_nan <- c(ind_nan.q, which(is.nan(gf_r))) # Remove NaN in global features
        gf<-cbind(gf_r,gf_q)
        if(length(ind_nan)>0){ gf<-gf[-ind_nan,] }
        if( (nrow(gf)>1) & (ncol(gf)==2)){ cos_value <- cosine(gf)[1,2] }
        rm(gf);rm(gf_r);rm(ind_nan);
      }  
      
      #save window start index according cos criteria and cos value
      if(cos_value>=threshold_cos_temp){
        pos_i <- c(pos_i, i)
        Cosine_threshold <- c(Cosine_threshold, cos_value)
      }
      # next window start 
      i <- i+step_threshold
    }
    threshold_cos_temp <- threshold_cos_temp-0.01
  }
  if(length(pos_i)==0){warning("No similar window looks appropriate for imputation")}
  return(pos_i)
}


#' @title Finding similar windows to a query
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description This function finds similar windows to a query consisting of a univariate vector.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.Finding_similar_window_univariate <- function(query, database, selected_qs, ...){
  id_similar=NULL;
  if (length(selected_qs)==0) return(id_similar)
  
  # Initialization
  T_gap <- length(query)
  Cost_dist <- c()
  pos_i <- selected_qs
  
  for (i in pos_i){
    k <- i+T_gap-1
    ref <- database[i:k]
    align <- dtw(query, ref, keep.internals=TRUE)
    cost <- align$normalizedDistance
    Cost_dist <- c(Cost_dist, cost)
  }
  
  min_cost <- min(Cost_dist)
  id_similar <- pos_i[which(Cost_dist==min_cost)]
  
  return(list(id_similar))
}



imputeDTWBI_QF<-function (data, t_gap, T_gap, S_Query, DTW_method = "DTW", threshold_cos = NULL, 
                          step_threshold = NULL, thresh_cos_stop = 0.8, verbose=F,...){
  
  method_options <- c("DTW", "DDTW", "AFBDTW", "dtw", "ddtw","afbdtw")
  if (DTW_method %in% method_options){
    if(verbose ==T) print(DTW_method)
  }else{ stop("Invalid DTW method")}
  
  IdGap <- t_gap:(t_gap + T_gap - 1)
  
  #if (sum(which(is.na(data[-IdGap]))) > 0) {
  #    stop("Dataset contains remaining NA outside main gap")}
  
  N <- length(data)
  if (S_Query >= 0.25*N) stop("Gap is too large to compute an appropriate imputation proposal")
  
  if (is.null(threshold_cos)){
    if (N > 10000) { threshold_cos = 0.9995
    }else { threshold_cos = 0.995}
  }
  
  if (is.null(step_threshold)){
    if (N > 10000) {step_threshold = 50
    }else if (N > 1000) {
      step_threshold = 10
    }else{ step_threshold = 2}
  }
  
  flag_tested=F;
  selected_qs <- NULL;
  id_similar_window=NULL;
  
  #test gap position : first half or second half in the data.
  #case 1: query after gap if t_gap<N/2
  if (t_gap < N/2) {
    flag_tested=T;
    query_a <- c()
    data_a <- data
    pos_start <- t_gap + T_gap
    ind <- pos_start:(pos_start + S_Query - 1)
    query_a <- data[ind]
    if (DTW_method == "DDTW" || DTW_method == "ddtw") {
      query_a <- local.derivative.ddtw(query_a)
      data_a <- local.derivative.ddtw(data_a)
      data_a[t_gap - 1] <- 0
      data_a[t_gap + T_gap] <- 0
    }
    i_start <- pos_start + max(S_Query,T_gap)
    i_finish <- length(data) - S_Query + 1
    
    selected_qs <- .DTW_threshold_global_univariate(query = query_a, 
                                                    database = data_a, i_start, i_finish, T_gap, T, step_threshold, 
                                                    threshold_cos, thresh_cos_stop, ...)
    if(length(selected_qs)>0){
      if (DTW_method == "AFBDTW" || DTW_method == "afbdtw") {
        id_similar_window <- .Finding_similar_window_univariate_AFBDTW(query = query_a, 
                                                                       database = data_a, selected_qs = selected_qs, 
                                                                       ...)
      }else {
        id_similar_window <- .Finding_similar_window_univariate(query = query_a, 
                                                                database = data_a, selected_qs = selected_qs, 
                                                                ...)
      }
      if(length(id_similar_window)>0){
        id_simwin_begin <- id_similar_window[[1]]#id du d?but de la fen?tre trouv?e
        id_simwin_end <- id_simwin_begin + S_Query - 1#id fin de la fen?tre trouv?e (donc S_Query)
        
        similar_query_dtw <- data[id_simwin_begin:id_simwin_end]#fen?tre
        id_imp_end <- id_simwin_begin - 1 #fin du remplacement
        id_imp_begin <- id_imp_end - T_gap + 1 #d?but du remplacement
        imp_value_dtw <- data[id_imp_begin:id_imp_end]#valeur ? prendre
        data_imputed_proposal = data
        data_imputed_proposal[t_gap:(t_gap + T_gap - 1)] <- imp_value_dtw
        imputation <- list(output_vector = data_imputed_proposal, 
                           input_vector = data, query = data[ind], pos_query = c(pos_start, 
                                                                                 (pos_start + T_gap - 1)), sim_window = similar_query_dtw, 
                           pos_sim_window = c(id_simwin_begin, id_simwin_end), 
                           imputation_window = imp_value_dtw, pos_imp_window = c(id_imp_begin, 
                                                                                 id_imp_end))
      }
    }
  }
  #case 2 and 3: query before if t_gap<N/2 or no selected qs for case 1
  if (t_gap >= N/2 | length(selected_qs)==0) {
    query_b <- c()
    pos_start <- t_gap - 1
    ind <- (pos_start - S_Query + 1):(pos_start)
    Researchbase_b <- data[1:pos_start]
    query_b <- data[ind]
    if (DTW_method == "DDTW" || DTW_method == "ddtw") {
      query_b <- local.derivative.ddtw(query_b)
      Researchbase_b <- local.derivative.ddtw(Researchbase_b)
    }
    i_start <- 1
    i_finish <- length(Researchbase_b) - 2 * max(S_Query,T_gap)
    selected_qs <- .DTW_threshold_global_univariate(query = query_b, 
                                                    database = Researchbase_b, i_start = i_start, i_finish = i_finish, T_gap,
                                                    F, step_threshold, threshold_cos, thresh_cos_stop)
    if(length(selected_qs)>0){
      if (DTW_method == "AFBDTW" || DTW_method == "afbdtw") {
        id_similar_window <- .Finding_similar_window_univariate_AFBDTW(query = query_b, 
                                                                       database = Researchbase_b, selected_qs = selected_qs)
      }else {
        id_similar_window <- .Finding_similar_window_univariate(query = query_b, 
                                                                database = Researchbase_b, selected_qs = selected_qs)
      }
      if(length(id_similar_window)>0){
        id_simwin_begin <- id_similar_window[[1]]
        id_simwin_end <- id_simwin_begin + S_Query - 1
        
        similar_query_dtw <- Researchbase_b[id_simwin_begin:id_simwin_end]
        id_imp_begin <- id_simwin_end + 1
        id_imp_end <- id_imp_begin + T_gap - 1
        imp_value_dtw <- data[id_imp_begin:id_imp_end]
        data_imputed_proposal = data
        data_imputed_proposal[t_gap:(t_gap + T_gap - 1)] <- imp_value_dtw
        imputation <- list(output_vector = data_imputed_proposal, 
                           input_vector = data, query = data[ind], pos_query = c((pos_start - 
                                                                                    T_gap + 1), pos_start), sim_window = similar_query_dtw, 
                           pos_sim_window = c(id_simwin_begin, id_simwin_end), 
                           imputation_window = imp_value_dtw, pos_imp_window = c(id_imp_begin, 
                                                                                 id_imp_end))
      }
    }
  }
  #case 4 
  if(flag_tested==F & length(selected_qs)==0){
    query_a <- c()
    data_a <- data
    pos_start <- t_gap + T_gap
    ind <- pos_start:(pos_start + S_Query - 1)
    query_a <- data[ind]
    if (DTW_method == "DDTW" || DTW_method == "ddtw") {
      query_a <- local.derivative.ddtw(query_a)
      data_a <- local.derivative.ddtw(data_a)
      data_a[t_gap - 1] <- 0
      data_a[t_gap + T_gap] <- 0
    }
    i_start <- pos_start + max(S_Query,T_gap)
    i_finish <- length(data) - S_Query + 1
    
    selected_qs <- .DTW_threshold_global_univariate(query = query_a, 
                                                    database = data_a, i_start, i_finish, T_gap, T, step_threshold, 
                                                    threshold_cos, thresh_cos_stop, ...)
    if(length(selected_qs)>0){
      if (DTW_method == "AFBDTW" || DTW_method == "afbdtw") {
        id_similar_window <- .Finding_similar_window_univariate_AFBDTW(query = query_a, 
                                                                       database = data_a, selected_qs = selected_qs, ...)
      }else{
        id_similar_window <- .Finding_similar_window_univariate(query = query_a, 
                                                                database = data_a, selected_qs = selected_qs, ...)
      }
      if(length(id_similar_window)>0){
        id_simwin_begin <- id_similar_window[[1]]
        id_simwin_end <- id_simwin_begin + S_Query - 1
        similar_query_dtw <- data[id_simwin_begin:id_simwin_end]
        id_imp_end <- id_simwin_begin - 1
        id_imp_begin <- id_imp_end - T_gap + 1
        imp_value_dtw <- data[id_imp_begin:id_imp_end]
        data_imputed_proposal = data
        data_imputed_proposal[t_gap:(t_gap + T_gap - 1)] <- imp_value_dtw
        imputation <- list(output_vector = data_imputed_proposal, 
                           input_vector = data, query = data[ind], pos_query = c(pos_start, (pos_start + T_gap - 1)), 
                           sim_window = similar_query_dtw, 
                           pos_sim_window = c(id_simwin_begin, id_simwin_end), 
                           imputation_window = imp_value_dtw, 
                           pos_imp_window = c(id_imp_begin, id_imp_end))
      }
    }
  }
  if(length(selected_qs)==0 | length(id_similar_window)==0){
    if(verbose) warning("no available window according global features or DTW cost")
    sigNA=rep(NA,T_gap)
    imputation <- list(output_vector = sigNA, 
                       input_vector = data, query = data[ind], 
                       pos_query = c(pos_start, (pos_start + T_gap - 1)), 
                       sim_window = sigNA, 
                       pos_sim_window = c(NA, NA), 
                       imputation_window = sigNA, 
                       pos_imp_window = c(NA, NA))
  }
  return(imputation)
}

#' calculate the trend of the series by smoothing it with a moving average
#' 
#' '
#' @param serie time series to smooth
#' @param rec number of time the smoothing is done
#' @param mov_win size of the window which calculates the moving average
#' @return trend of the serie
lissage <- function(serie, rec=1, mov_win=NULL){
  if(is.null(mov_win)){mov_win=floor(length(serie)/20)}
  smooth=rep(0,length(serie))
  out=serie
  for (k in 1:rec){
    smooth[1:mov_win]=NA
    smooth[(length(serie)-mov_win):length(serie)]=NA
    
    beg=mov_win+1
    end=length(serie)-beg
    
    for(i in 1:(2*mov_win+1)){
      smooth[beg:end]=smooth[beg:end]+out[i:(length(serie)-2*mov_win-2+i)]/(2*mov_win+1)
    }
    smooth[1]=serie[1]
    smooth[length(serie)]=serie[length(serie)]
    smooth=na.interp(smooth)[1:length(serie)]
    
    out=smooth
    smooth=rep(0,length(serie))
  }
  out[1:mov_win]=NA
  out[(length(serie)-mov_win):length(serie)]=NA
  return(out)
}

#' calculate the trend of the series by smoothing it with a moving average.
#' a small number of missing data is ignored when calculating the trend
#' 
#' '
#' version 18.02.2023
#' @param serie time series to smooth
#' @param rec number of time the smoothing is done
#' @param mov_win size of the window which calculates the moving average
#' @return trend of the serie
lissage2 <- function(serie, rec=1, mov_win=NULL){
  nbObs=length(serie)
  if(is.null(mov_win)){mov_win=floor(nbObs/20)}
  smooth=rep(0,nbObs)
  out=serie
  for (k in 1:rec){
    if(sum(is.na(out))<nbObs){
      smooth[1:mov_win]=NA
      smooth[(nbObs-mov_win):nbObs]=NA
    
      beg=mov_win+1
      end=nbObs-beg
    
      for(i in beg:end){
        gap=((i-mov_win):(i+mov_win))
        if (sum(is.na(out[gap]))<=50){
          smooth[i]=mean(na.omit(out[gap]))
        }
        else{smooth[i]=NA}
      }
    
      if(sum(is.na(smooth))!=nbObs){
        smooth=na.interp(smooth)[1:nbObs]
      }
      out=smooth
      smooth=rep(0,nbObs)
    }
  }
  out[1:mov_win]=NA
  out[(nbObs-mov_win):nbObs]=NA
  return(out)
}

#' decompose the time series into 3 sub-series : the trend, the cycle and the residue
#' 
#' '
#' @author Marco Hanocq, Emilie Poisson Caillault
#' @version 18.02.2023
#' @param serie time series to decompose
#' @param mov_win size of the window which calculates the moving average for the trend
#' @param .seasonak if known, what is the period of the time series
#' @param mode select either additive or multiplicative decomposition
#' @param ratio between 0 and 1 to accept to compute trend/cycle. 
#' @return a list which countains the original series and the 3 sub-series mentioned above
compute.decompose <- function(serie,mov_win=NULL,.seasonal=NULL, mode = "additif",ratio=0.3){
  nbObs=length(serie)
  nbNA=sum(is.na(serie))
  if(is.null(mov_win)){mov_win=floor(nbObs/20)}
  res=NULL
  if(mode=="additif"){
    res=list("x"=serie,"trend"=rep(0,nbObs),"seasonal"=rep(0,nbObs),"random"=serie)
  }else{
    res=list("x"=serie,"trend"=rep(1,nbObs),"seasonal"=rep(1,nbObs),"random"=serie)
  }
  
  if(nbNA/nbObs>ratio){
    trendS=lissage2(serie,2,mov_win)
    if(sum(is.na(trendS))<nbObs){
      res$trend=trendS
    }
     rm(trendS)
    
     if (mode=="additif") {
        temp=serie-res$trend
     }else {temp= serie/res$trend}
  
     seasonal=res$seasonal
    if(is.null(.seasonal)){
      seasonal=compute.bestperiod(serie,temp)
    }else{
      for(i in 1:.seasonal){
        seasonal[seq(i,nbObs,.seasonal)]=mean(na.omit(temp[seq(i,nbObs,.seasonal)]))
      }
    }
  
    if (mean(seasonal)!=0){res$seasonal=seasonal; rem(seasonal)}  
    if (mode=="additif") {
      res$random=res$x-res$trend-res$seasonal
    }else {res$random=res$x/res$trend/res$seasonal}
  }
  return(res)
}

#' compute the period of time series without trend
#' 
#' '
#' @param serie original series
#' @param untrend original series whithout its trend
#' @return a list which countains the original series and the 3 sub-series mentioned above
compute.bestperiod <- function(serie,untrend){
  nbObs=length(serie)
  temp1=rep(0,nbObs)
  m=0
  p=1
  for(period in 12:nbObs/10){
    for(i in 1:period){
      temp1[seq(i,nbObs,period)]=mean(na.omit(untrend[seq(i,nbObs,period)]))
    }
    if(sum(is.na(temp1))!=0){
      print(is.na(temp1))
      print(i)
    }
    if(max(temp1)>m){
      m=max(temp1)
      res=temp1
      p=period
    }
  }
  #print(period)
  return(res)
}

#' complete a time series by decomposing it in sub-series and completing the sub-series
#' 
#' '
#' @author  Marco Hanocq, Emilie Poisson Caillault
#' @version 10.02.2023 
#' @param serie time series to complete
#' @param acceptedHole size of the holes that are considered too big to complete
#' @param smallHole size of the holes that are completed by adaptative moving average
#' @param S_Query size of the query
#' @param .seasonal if the periodic component of the time serie is known, which is it
#' @param mov_win size of the window use in the moving average when computing the trend
#' @param mode additive decomposition or multiplicative decomposition
#' @return completed time serie
#' 
#' @example : compute.decompose_holed_serie(imputeTS::tsNH4,500,12,12*40)
compute_decompose_holed_series<-function(serie, acceptedHole, smallHole, S_Query, .seasonal = NULL, mov_win = NULL, mode = "additif"){
  nbObs=length(serie)
  if(sum(is.na(serie))==0){return(serie)}
  if(is.null(mov_win)){mov_win=floor(nbObs/20)}
  
  dd=completion.singleHole(serie)
  s1=completion.adaptativeMovingAverage(dd,acceptedHole = smallHole,mult = 10,confidence = 1,maxWindow = 10 ,V = F);
  s2=completion.adaptativeMovingAverage(dd[nbObs:1],acceptedHole = smallHole,mult = 10,confidence = 1,maxWindow = 10,V = F);
  s2=s2[nbObs:1]
  dd=(s1+s2)/2;
  rm(s1);rm(s2)
  
  temp1=dd[1:mov_win]
  temp2=dd[(nbObs-mov_win):nbObs]
  
  dd_res=compute.decompose(dd, mov_win = mov_win, .seasonal = .seasonal, mode)
  
  dd$random=compute.DTWBI_QF(dd$random,acceptedHole,smallHole,S_Query,verbose = F)
  dd$random=dd$random[1:length(dd$random)]
  
  if(mode == "additif"){
    dd$x=dd$trend+dd$seasonal+dd$random
  }else{
    dd$x=dd$trend*dd$seasonal*dd$random
  }
  dd$x[1:mov_win]=temp1
  dd$x[(nbObs-mov_win):nbObs]=temp2
  
  return(dd)
}






