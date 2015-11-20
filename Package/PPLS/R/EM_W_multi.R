#' PPLS: Probabilistic Partial Least Squares
#'
#' This package implements the Partial Least Squares method in a probabilistic framework.
#' @author
#' Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}),
#' Jeanine Houwing-Duistermaat (\email{J.J.Houwing@@lumc.nl}),
#' Geurt Jongbloed (\email{G.Jongbloed@@tudelft.nl}),
#' Szymon Kielbasa (\email{S.M.Kielbasa@@lumc.nl}),
#' Hae-Won Uh (\email{H.Uh@@lumc.nl}).
#'
#' Maintainer: Said el Bouhaddani (\email{s.el_bouhaddani@@lumc.nl}).
#'
#' @section Usage:
#' Just do PPLS(Datamatrix_1,Datamatrix_2,nr_of_components) to have a quick and dirty fit, modify the default arguments if necessary.
#'
#' @docType package
#' @name PPLS
#' @keywords Probabilistic-PLS
#' @import Rcpp
#' @import ggplot2
#' @import RcppEigen
#' @import O2PLS
NULL

#' Performs one EM step
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param W. Numeric matrix.
#' @param C. Numeric matrix.
#' @param B_T. Numeric
#' @param sigX. Numeric
#' @param sigY. Numeric
#' @param sigH. Numeric
#' @param sigT. Numeric
#' @return A list with updated values for\itemize{
#'    \item{W}{Matrix}
#'    \item{C}{Matrix}
#'    \item{B}{Matrix}
#'    \item{sighat}{Vector containing updated sigX and sigY}
#'    \item{siglathat}{Vector containing sigH and sigT}
#'  }
#' @details This function passes its arguments to EMstepC (a C++ function, see \link{Rcpp}), which returns the expected sufficient statistics.
#' The maximization is done afterwards in this function. This may become a full C++ function later.
#'
#' @export
EMstep_W <- function(X , Y , W.=W, C.=C, B_T.=B_T,
                     sigX.=sigX,sigY.=sigY,sigH.=sigH,sigT.=sigT)
{
  W. = as.vector(W.)
  C. = as.vector(C.)
  B_T. = B_T.[1]
  sigX. = sigX.[1]
  sigY. = sigY.[1]
  sigH. = sigH.[1]
  sigT. = sigT.[1]

  g = sigT.^2*B_T.^2 + sigH.^2
  Kw = sigT.^2 - sigT.^4*B_T.^2/sigY.^2 + sigT.^4*B_T.^2*g/(sigY.^2*(g+sigY.^2))
  Kc = g - sigT.^4*B_T.^2/sigX.^2 + sigT.^6*B_T.^2/(sigX.^2*(sigT.^2+sigX.^2))
  Kwc = sigT.^2*B_T./(sigX.^2*sigY.^2) - Kc*sigT.^2*B_T./(sigX.^2*sigY.^2*(Kc+sigY.^2)) -
    sigT.^4*B_T./(sigX.^2*sigY.^2*(sigT.^2+sigX.^2)) +
    Kc*sigT.^4*B_T./(sigX.^2*sigY.^2*(Kc+sigY.^2)*(sigT.^2+sigX.^2))
  c1 = Kw / (sigX.^2*(Kw + sigX.^2));
  c3 = Kc / (sigY.^2*(Kc + sigY.^2));
  c2 = Kwc;

  return(EMstepC_fast(W.,C.,B_T.,X,Y,sigX.,sigY.,sigH.,sigT.,c1,c2,c3))

}


#' Performs PPLS fit in one direction.
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param EMsteps strictly positive integer. Denotes the maximum number of EM steps to make.
#' @param atol Double, convergence criterium for the log-likelihood
#' @param initialGuess A string. Choose "o2m", "random" or "equal" depending on what type of initial guess you want.
#' @param critfunc Function measuring the increment value, i.e. f(L_{i+1} - L_{i}). Usually f == I or f == abs.
#' @param print_logvalue Logical. Should we print log-likelihood values for each step?
#' @return A list with estimates for\itemize{
#'    \item{W}{Matrix}
#'    \item{C}{Matrix}
#'    \item{B}{Matrix}
#'    \item{sig}{Vector containing updated sigX, sigY, sigH and sigT}
#'    \item{logvalue}{Vector with last (or all, if print_logvalue==TRUE) log-likelihood(s)}
#'    \item{Last_increment}{Double, value of last increment.}
#'    \item{Number_steps}{Integer, number of steps needed to stop.}
#'  }
#' @details This function estimates loadings and variances for one direction.
#'
#' @export
PPLSi <- function(X,Y,EMsteps=1e2,atol=1e-4,initialGuess=c("o2m","random","equal"),critfunc=function(x){x},print_logvalue=FALSE)
{
  a = 1
  p = ncol(X)
  q = ncol(Y)
  N = nrow(X)
  initialGuess = match.arg(initialGuess)
  if(initialGuess == "o2m"){
    sim=o2m(X,Y,1,0,0)
    Wnw=sim$W.;Cnw=sim$C.;Bnw=sim$B_T.[1];
    signw=sqrt(c(ssq(sim$E)/N/p,ssq(sim$F)/N/q));
    siglatnw=sqrt(c(ssq(sim$H_U),ssq(sim$Tt))/N)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  } else if(initialGuess == "random"){
    Wnw=orth(runif(p));Cnw=orth(runif(q)); Bnw = rchisq(1,1); siglatnw =  rchisq(2,100)/100; signw = rchisq(2,10)/100
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }else if(initialGuess == "equal"){
    Wnw=orth(rep(1,p));Cnw=orth(rep(1,q)); Bnw = 1; siglatnw =  c(1,1); signw = c(1/p,1/q)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }

  logvalue=1:(EMsteps+1)*NA

  logvalue[1]=logl_W(X,Y,Wnw,Cnw,Bnw,signw[1],signw[2],siglatnw[1],siglatnw[2])

  for(i in 1:EMsteps){
    if(any(signw < 100*.Machine$double.eps)){
      return(list(W=NA,C=NA,B=NA,sig=NA,logvalue=NA,Last_increment = NA, Number_steps = i))
    }
    fit=EMstep_W(X,Y,Wnw,Cnw,Bnw,signw[1],signw[2],siglatnw[1],siglatnw[2])
    Bnw = (fit$B)
    Wnw = (fit$W)
    #PYonw = (fit$PY)
    Cnw = (fit$C)
    #PXonw = (fit$PX)
    signw = fit$sigh
    siglatnw = fit$sigl
    logvalue[i+1]=logl_W(X,Y,Wnw,Cnw,Bnw,signw[1],signw[2],siglatnw[1],siglatnw[2])
    if(critfunc(logvalue[i+1]-logvalue[i])<atol){break}
    #on.exit(return(list( W=c(Wnw),C=c(Cnw),B=Bnw,sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = logvalue[i]-logvalue[i-1], Number_steps = i)))
  }
  last_incr = logvalue[i+1]-logvalue[i]
  if(any(diff(logvalue[0:i+1])<0)){warning("Not monotone")}
  if(critfunc(logvalue[i+1]-logvalue[i])>=atol){warning(paste("Not converged, last increment was",critfunc(logvalue[i+1]-logvalue[i])))}
  #if(!print_logvalue){logvalue=NA}
  ret2 = list( W=c(Wnw),C=c(Cnw),B=c(Bnw),sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = last_incr, Number_steps = i)
  return(ret2)
}


#' Performs PPLS fit in sequentially multiple directions.
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param nr_comp Strictly positive integer. The number of PPLS components to fit.
#' @param EMsteps strictly positive integer. Denotes the maximum number of EM steps to make.
#' @param atol Double, convergence criterium for the log-likelihood
#' @param initialGuess A string. Choose "o2m", "random" or "equal" depending on what type of initial guess you want.
#' @param critfunc Function measuring the increment value, i.e. f(L_{i+1} - L_{i}). Usually f == I or f == abs.
#' @param print_logvalue Logical. Should we print log-likelihood values for each step?
#' @return A list with estimates for\itemize{
#'    \item{W}{Matrix, each column is a fit.}
#'    \item{C}{Matrix, each column is a fit.}
#'    \item{B}{Matrix, each column is a fit.}
#'    \item{sig}{Matrix containing updated sigX, sigY, sigH and sigT. Columns are named.}
#'    \item{Other_output}{\itemize{
#'          \item{Last_increment}{.}
#'          \item{Number_steps}{.}
#'          \item{Loglikelihoods}{.}
#'          }
#'    }
#'  }
#' @details This function estimates loadings and variances for multiple orthogonal directions.
#'
#' @export
PPLS <- function(X,Y,nr_comp=1,EMsteps=1e2,atol=1e-4,initialGuess=c("equal","o2m","random"),critfunc=function(x){x},print_logvalue=FALSE)
{
  initialGuess = match.arg(initialGuess)
  a = nr_comp
  p = ncol(X)
  q = ncol(Y)
  N = nrow(X)
  if(p+q>5000 && initialGuess == 'o2m'){initialGuess = 'equal';warning("p or q too large for 'o2m' initial guess, we chose 'equal' instead.")}
  Wn = matrix(NA,p,a)
  Cn = matrix(NA,q,a)
  Bn = sigXn = sigYn = sigHn = sigTn = 1:a*NA
  other_output = list()
  Xc = X
  Yc = Y
  for(i in 1:a){
    fit = PPLSi(Xc,Yc,EMsteps,atol,initialGuess,critfunc,print_logvalue=TRUE);
    if(is.na(fit$B[1])){
      warning(paste("From component",i,"on the residuals are of rank <",1e-14,"and calculations are stopped."))
      i = i-1;
      Wn = Wn[,1:i]; Cn = Cn[,1:i]; Bn = Bn[1:i]; sigXn = sigXn[1:i]; sigYn = sigYn[1:i]; sigHn = sigHn[1:i]; sigTn = sigTn[1:i];
      break
    }
    Wn[,i] = c(fit$W);
    Cn[,i] = c(fit$C);
    Bn[i] = c(fit$B);
    sigXn[i] = fit$sig[1]
    sigYn[i] = fit$sig[2]
    sigHn[i] = fit$sig[3]
    sigTn[i] = fit$sig[4]
    Xc = Xc - (Xc%*%fit$W)%*%t(fit$W)
    Yc = Yc - (Yc%*%fit$C)%*%t(fit$C)
    other_output$Last_increment[i] = fit$Last
    other_output$Number_steps[i] = fit$Number
    other_output$Loglikelihoods[i] = logl_W(X,Y,Wn[,1:i],Cn[,1:i],diag(Bn,i),sigXn[i],sigYn[i],sigHn[i],diag(sigTn[1:i],i))
  }
  ret2 = list( W=Wn,C=Cn,B=Bn,sig=cbind(sigX=sigXn,sigY=sigYn,sigH=sigHn,sigT=sigTn),Other_output=other_output)
  class(ret2) <- "PPLS"
  return(ret2)
}

#' Calculates loglikelihood for PPLS model.
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param W. Numeric matrix.
#' @param C. Numeric matrix.
#' @param B_T. Numeric matrix.
#' @param sigX. Numeric matrix.
#' @param sigY. Numeric matrix.
#' @param sigH. Numeric matrix.
#' @param sigT. Numeric matrix.
#' @return The log-likelihood (i.e. scalar) of the parameters given the observed data.
#'
#' @details This function passes its arguments to LoglC_fast (a C++ function, see \link{Rcpp}).
#'
#' @export
logl_W <- function(X , Y , W.=W, C.=C, B_T.=B_T,
                   sigX.=sigX,sigY.=sigY,sigH.=sigH,sigT.=sigT)
{
  W. = as.matrix(W.)
  C. = as.matrix(C.)
  p = nrow(W.)
  q = nrow(C.)
  a = ncol(W.)

  B_T. = t(B_T.)
  sigX. = sigX.[1]
  sigY. = sigY.[1]
  sigH. = sigH.[1]
  sigT. = t(sigT.)

  g = diag(sapply(1:a,function(i)sqrt(sigT.[i,i]^2%*%B_T.[i,i]^2 + sigH.^2)),a)
  Kw = diag(sapply(1:a,function(i)sigT.[i,i]^2 - sigT.[i,i]^4*B_T.[i,i]^2/sigY.^2 + sigT.[i,i]^4*B_T.[i,i]^2*g[i,i]^2/(sigY.^2*(g[i,i]^2+sigY.^2))),a)
  Kc = diag(sapply(1:a,function(i)g[i,i]^2 - sigT.[i,i]^4*B_T.[i,i]^2/sigX.^2 + sigT.[i,i]^6*B_T.[i,i]^2/(sigX.^2*(sigT.[i,i]^2+sigX.^2))),a)
  Kwc = diag(sapply(1:a,function(i) sigT.[i,i]^2*B_T.[i,i]/(sigX.^2*sigY.^2) - Kc[i,i]*sigT.[i,i]^2*B_T.[i,i]/(sigX.^2*sigY.^2*(Kc[i,i]+sigY.^2)) -
                      sigT.[i,i]^4*B_T.[i,i]/(sigX.^2*sigY.^2*(sigT.[i,i]^2+sigX.^2)) +
                      Kc[i,i]*sigT.[i,i]^4*B_T.[i,i]/(sigX.^2*sigY.^2*(Kc[i,i]+sigY.^2)*(sigT.[i,i]^2+sigX.^2)) ),a)
  c1 = diag(sapply(1:a,function(i) Kw[i,i] / (sigX.^2*(Kw[i,i] + sigX.^2)) ),a)
  c3 = diag(sapply(1:a,function(i) Kc[i,i] / (sigY.^2*(Kc[i,i] + sigY.^2)) ),a)
  c2 = Kwc

  return(loglC_fast(W.,C.,X,Y,sigX.,sigY.,diag(sigT.)^2,diag(c1),diag(c2),diag(c3),diag(Kc)))
}

#' Print function for class PPLS.
#'
#' @param x A PPLS fit (an object of class PPLS)
#' @param perc Display relative percentages yes/no?
#' @param ... For consistency
#' @return Invisible, the outputted table
#'
#' @details This function shows the absolute/relative variances of each component.
#'
#' @export
print.PPLS <- function (x,perc=FALSE,...)
{
  p = nrow(x$W)
  q = nrow(x$C)
  outp = sapply(1:nrow(x$sig), function(i) {
    with(x, {
      c(i,sum(sig[1:i, 4]^2)/(perc*(sum(sig[1:i, 4]^2) + p * x$sig[i, 1]^2) + (1 - perc)),
        sum(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2)/(perc * (sum(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2) + q * sig[i, 2]^2) +(1 - perc)),
        sig[i, 3]^2 / (perc*(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2) + (1-perc)),
        c(0, diff(x$Oth$Log))[i],x$Oth$Nu[i], signif(x$Oth$Last[i], 3))
    })
  })
  outp = as.data.frame(t(outp))
  names(outp) = c("LV", ifelse(perc, "ssq(T)/ssq(X)", "ssq(T)"),
                  ifelse(perc, "ssq(U)/ssq(Y)", "ssq(U)"), "sigma_H", "log LR",
                  "#steps", "last incr")
  print(outp)
  return(invisible(outp))
}

#' Plot function for class PPLS.
#'
#' @param x A PPLS fit (an object of class PPLS)
#' @param XorY Plot loadings for X or Y?
#' @param i Positive integer. The i-th loading/score will be plotted on the first axis
#' @param j Positive integer. The j-th loading/score will be plotted on the second axis
#' @param use_ggplot2 Logical. Use ggplot2?
#' @param ... For consistency
#' @return If ggplot2 is active, returns the ggplot object. Else NULL.
#'
#' @details This function plots two loadings/scores.
#'
#' @export
plot.PPLS <- function (x, XorY = c("X", "Y"), i = 1, j = NULL, use_ggplot2=TRUE,...)
{
  fit = list()
  XorY = match.arg(XorY)
  if(XorY == "X"){fit$load = x$W[,c(i,j)]}else{fit$load = x$C[,c(i,j)]}
  p = nrow(as.matrix(fit$load))
  if(is.null(j)){
    fit$load = cbind(1:p,fit$load)
    colnames(fit$load) = c("index",paste("loadings",i))
  }else{
    colnames(fit$load) = c(paste("loadings",i),paste("loadings",j))
  }
  if (use_ggplot2) {
    plt = with(fit, qplot(x = load[, 1], y = load[, 2], label = 1:p,
                          geom = "text", xlab = colnames(load)[1], ylab = colnames(load)[2]))
    plt = plt + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
    print(plt)
    return(plt)
  }
  else {
    with(fit, {
      plot(load[, 1], load[, 2], type = "n")
      text(load[, 1], load[, 2])
    })
    abline(v=0,h=0)
    #lines((2/sqrt(p) + 2/p) * cos(seq(0, 2 * pi, length.out = 1000)),
    #      (2/sqrt(p) + 2/p) * sin(seq(0, 2 * pi, length.out = 1000)))
  }
}

#' Get scores from PPLS fit.
#'
#' @param fit A PPLS fit (an object of class PPLS)
#' @param X Matrix
#' @param Y Matrix
#' @param subset vector of positive integers denoting with components you want.
#' @return Both the X and Y scores concatenated.
#'
#' @details This function plots two loadings/scores.
#'
#' @export
scores.PPLS <- function (fit, X, Y, subset = NULL)
{
  if (is.null(subset)) {
    subset = 1:ncol(fit$W)
  }
  if (length(subset) == 1) {
    return(c(X %*% fit$W[, subset], Y %*% fit$C[, subset]))
  }
  return(rbind(X %*% fit$W[, subset], Y %*% fit$C[, subset]))
}
