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
#' @name Probabilistic-PLS
#' @keywords Probabilistic-PLS
#' @import Rcpp
#' @import ggplot2
#' @import RcppEigen
#' @import O2PLS
#' @import magrittr
#' @useDynLib PPLS
NULL

#' Performs one EM step (use \link{PPLS} for fitting a PPLS model)
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
#'    \item{W }{Matrix}
#'    \item{C }{Matrix}
#'    \item{B }{Matrix}
#'    \item{sighat }{Vector containing updated sigX and sigY}
#'    \item{siglathat }{Vector containing sigH and sigT}
#'  }
#' @details This function passes its arguments to EMstepC (a C++ function, see \link{Rcpp}), which returns the expected sufficient statistics.
#' The maximization is done afterwards in this function. This may become a full C++ function later.
#'
#' @export
EMstep_W <- function(X , Y , W., C., B_T.,
                     sigX.,sigY.,sigH.,sigT.)
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

#' Defines an appropriate list with chosen constraints
#'
#' @param constraints List. Should contain one or more numeric elements with name(s) from W, C, B, sigE, sigF, sigH and sigT. Defaults to \code{NULL} is no constraints are given.
#'
#' @return A named list with elements either the value of the constraint or \code{NULL}.
#' @examples
#' fconstraint(list(sigH=.1))
#' fconstraint(list(B=3,sigE=.5))
#' @export
fconstraint <- function(constraints=NULL){
  outp = list(W=NULL,C=NULL,B=NULL,sigE=NULL,sigF=NULL,sigH=NULL,sigT=NULL)
  if(is.null(constraints)){return(outp)}
  names_constraints = names(constraints)
  defined_constraints = which(names(outp) %in% names_constraints)
  outp[defined_constraints] = constraints[which(names_constraints %in% names(outp))]
  return(outp)
}

#' Performs PPLS fit in one direction. (use \link{PPLS} for fitting a PPLS model)
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param EMsteps strictly positive integer. Denotes the maximum number of EM steps to make.
#' @param atol Double, convergence criterium for the log-likelihood
#' @param initialGuess A string. Choose "o2m", "random" or "equal" depending on what type of initial guess you want.
#' @param customGuess A list with named components W,C,B,sigE,sigF,sigH,sigT.
#' @param critfunc Function measuring the increment value, i.e. f(L_{i+1} - L_{i}). Usually f == I or f == abs.
#' @param constraints List. Ideally being constructed with \link{fconstraint}.
#' @return A list with estimates for\itemize{
#'    \item{W}{Matrix}
#'    \item{C}{Matrix}
#'    \item{B}{Matrix}
#'    \item{sig}{ Vector containing updated sigX, sigY, sigH and sigT}
#'    \item{logvalue}{ Vector with last log-likelihood(s)}
#'    \item{Last_increment}{ Double, value of last increment.}
#'    \item{Number_steps}{ Integer, number of steps needed to stop.}
#'  }
#' @details This function estimates loadings and variances for one direction.
#'
#' @export
PPLSi <- function(X,Y,EMsteps=1e2,atol=1e-4,initialGuess=c("equal","o2m","random","custom"),customGuess=NULL,critfunc=function(x){x},
                  constraints=fconstraint())
{
  stopifnot(length(constraints) == 7,all(names(constraints) == c("W","C","B","sigE","sigF","sigH","sigT")))
  a = 1
  p = ncol(X)
  q = ncol(Y)
  N = nrow(X)
  initialGuess = match.arg(initialGuess)
  if(!is.null(customGuess)){initialGuess = "custom"}
  if(initialGuess == "o2m"){
    sim=o2m(X,Y,1,0,0)
    Wnw=sim$W.;Cnw=sim$C.;Bnw=sim$B_T.[1];
    signw=sqrt(c((ssq(X)-ssq(sim$Tt))/N/p,(ssq(Y)-ssq(sim$U))/N/q));
    siglatnw=sqrt(c(ssq(sim$U)-ssq(sim$Tt*Bnw),ssq(sim$Tt))/N)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  } else if(initialGuess == "random"){
    Wnw=orth(runif(p));Cnw=orth(runif(q)); Bnw = rchisq(1,1); siglatnw =  rchisq(2,100)/100; signw = rchisq(2,10)/100
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }else if(initialGuess == "equal"){
    Wnw=orth(rep(1,p));Cnw=orth(rep(1,q)); Bnw = 1; siglatnw =  c(1,1); signw = c(1/p,1/q)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }else if(initialGuess == "custom"){
    stopifnot('list' %in% class(customGuess))
    Wnw = customGuess$W;Cnw=customGuess$C;Bnw=customGuess$B;siglatnw=with(customGuess,c(sigH,sigT));signw=with(customGuess,c(sigE,sigF))
  }
  Wnw = with(constraints,if(is.numeric(W)) W else Wnw )
  Cnw = with(constraints,if(is.numeric(C)) C else Cnw )
  Bnw = with(constraints,if(is.numeric(B)) B else Bnw )
  signw = with(constraints,c(if(is.numeric(sigE)) sigE else signw[1] , if(is.numeric(sigF)) sigF else signw[2]))
  siglatnw = with(constraints,c(if(is.numeric(sigH)) sigH else siglatnw[1] , if(is.numeric(sigT)) sigT else siglatnw[2]))

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

    Wnw = with(constraints,if(is.numeric(W)) W else Wnw )
    Cnw = with(constraints,if(is.numeric(C)) C else Cnw )
    Bnw = with(constraints,if(is.numeric(B)) B else Bnw )
    signw = with(constraints,c(if(is.numeric(sigE)) sigE else signw[1] , if(is.numeric(sigF)) sigF else signw[2]))
    siglatnw = with(constraints,c(if(is.numeric(sigH)) sigH else siglatnw[1] , if(is.numeric(sigT)) sigT else siglatnw[2]))

    logvalue[i+1]=logl_W(X,Y,Wnw,Cnw,Bnw,signw[1],signw[2],siglatnw[1],siglatnw[2])
    if(critfunc(logvalue[i+1]-logvalue[i])<atol){break}
    #on.exit(return(list( W=c(Wnw),C=c(Cnw),B=Bnw,sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = logvalue[i]-logvalue[i-1], Number_steps = i)))
  }
  last_incr = logvalue[i+1]-logvalue[i]
  if(any(diff(logvalue[0:i+1])<0)){warning("Not monotone")}
  if(critfunc(logvalue[i+1]-logvalue[i])>=atol){warning(paste("Not converged, last increment was",critfunc(logvalue[i+1]-logvalue[i])))}
  ret2 = list( W=c(Wnw),C=c(Cnw),B=c(Bnw),sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = last_incr, Number_steps = i)
  return(ret2)
}


#' Performs PPLS fit.
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param nr_comp Strictly positive integer. The number of PPLS components to fit.
#' @param EMsteps strictly positive integer. Denotes the maximum number of EM steps to make.
#' @param atol Double, convergence criterium for the log-likelihood
#' @param initialGuess A string. Choose "o2m", "random" or "equal" depending on what type of initial guess you want.
#' @param customGuess A list with named components W,C,B,sigE,sigF,sigH,sigT.
#' @param critfunc Function measuring the increment value, i.e. f(L_{i+1} - L_{i}). Usually f == I or f == abs.
#' @param constraints List. Each element of \code{constraints} contains a named list, ideally being constructed with \link{fconstraint}. The numer of elements in \code{constraints} must equal nr_comp.
#' @return A list with estimates for\itemize{
#'    \item{W }{Matrix, each column is a fit.}
#'    \item{C }{Matrix, each column is a fit.}
#'    \item{B }{Matrix, each column is a fit.}
#'    \item{sig }{Matrix containing updated sigX, sigY, sigH and sigT. Columns are named.}
#'    \item{Other_output }{\itemize{
#'          \item{Last_increment}{ Last likelihood increment}
#'          \item{Number_steps}{ Number of EM steps}
#'          \item{Loglikelihoods}{ Likelihood at final step}
#'          }
#'    }
#'  }
#' @details This function estimates loadings and variances for multiple orthogonal directions. The loadings are W and C.
#' The X-scores can be retrieved via \code{X \%*\% fit$W} where fit is a PPLS fit, or with the function \code{\link{scores.PPLS}}.
#' The o2m initial guess fits a PLS model using the function \code{O2PLS::o2m} and uses the result as starting values.
#' If for some reason the likelihood is not monotone, check your input data or (re)try with \code{initialGuess = 'random'}.
#'
#' Use \code{plot(fit)} to plot the first two loadings. Use \code{print(fit,perc=TRUE)} for variances in percentages.
#'
#' @examples
#' exX = scale(matrix(rnorm(100*10),100,10))
#' exY = scale(matrix(rnorm(100*12),100,12))
#' PPLS(X = exX, Y = exY, nr_comp = 3, EMsteps = 1e4)
#' PPLS(X = exX, Y = exY, nr_comp = 3, EMsteps = 1e4, initialGuess = "random")
#'
#' # devtools::install_github("selbouhaddani/O2PLS")
#' if(require(O2PLS)){
#'  exinitGuess = list(W = orth(1:10), C = orth(1:12), B = 0.1,
#'                      sigE = 1,sigF = 1,sigH = 1,sigT = 0.1)
#'  PPLS(X = exX, Y = exY, nr_comp = 1, EMsteps = 1e4,
#'        initialGuess = "custom", customGuess = exinitGuess)
#' }
#' exconstraints = list(fconstraint(list(B = 1)) , fconstraint(list(sigT = 1)))
#' PPLS(X = exX, Y = exY, nr_comp = 2, EMsteps = 1e4, constraints = exconstraints)
#' @export
PPLS <- function(X,Y,nr_comp=1,EMsteps=1e2,atol=1e-4,initialGuess=c("equal","o2m","random","custom"),customGuess=NULL,critfunc=function(x){x},
                 constraints = lapply(1:nr_comp,fconstraint))
{
  if(!is.null(row.names(X))&!is.null(row.names(Y))){
    if(!all.equal(row.names(X),row.names(Y))){warning("specified rownames not equal!")}
  }
  X = as.matrix(X)
  Y = as.matrix(Y)
  stopifnot(nrow(X) == nrow(Y),ncol(X)>=nr_comp,ncol(Y)>=nr_comp)
  stopifnot(length(nr_comp)==1,length(EMsteps)==1,length(atol)==1)
  if(nr_comp<=0){stop("#components must be >0")}
  if(length(constraints) != nr_comp){stop("There should be a list of constraints for each component, see ?PPLS.")}
  initialGuess = match.arg(initialGuess)
  a = nr_comp
  p = ncol(X)
  q = ncol(Y)
  N = nrow(X)
  #if(p+q>5000 && initialGuess == 'o2m'){initialGuess = 'equal';warning("p or q too large for 'o2m' initial guess, we chose 'equal' instead.")}
  Wn = matrix(NA,p,a)
  Cn = matrix(NA,q,a)
  Bn = sigXn = sigYn = sigHn = sigTn = 1:a*NA
  other_output = list()
  Xc = X
  Yc = Y
  for(i in 1:a){
    if(any(sapply(constraints[[i]],is.numeric))){cat("constraints used on",names(constraints[[i]])[which(sapply(constraints[[i]],is.numeric))],"\n\n")}
    fit = PPLSi(X = Xc,Y = Yc,EMsteps = EMsteps,atol = atol,initialGuess = initialGuess,
                customGuess = customGuess,critfunc = critfunc,constraints = constraints[[i]]);
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
logl_W <- function(X , Y , W., C., B_T.,
                   sigX.,sigY.,sigH.,sigT.)
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
#' @param digits How many digits should be retained?
#' @param ... For consistency
#' @return Invisible data.frame, the outputted table
#'
#' @details This function shows the absolute/relative variances of each component.
#'
#' @export
print.PPLS <- function (x,perc=TRUE,digits=3,...)
{
  p = nrow(x$W)
  q = nrow(x$C)
  outp = sapply(1:nrow(x$sig), function(i) {
    with(x, {
      c(i,sum(sig[1:i, 4]^2)/(perc*(sum(sig[1:i, 4]^2) + p * x$sig[i, 1]^2) + (1 - perc)), #sigT
        sum(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2)/(perc * (sum(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2) + q * sig[i, 2]^2) +(1 - perc)), #sigU
        sig[i, 3]^2 / (perc*sum(sig[1:i,4]^2 * B[1:i]^2 + sig[i, 3]^2) + (1-perc)), #sigH
        c(0, diff(x$Oth$Log))[i],x$Oth$Nu[i], signif(x$Oth$Last[i], 3)) #the rest
    })
  })
  outp = as.data.frame(t(outp))
  names(outp) = c("LV", ifelse(perc, "ssq(T)/ssq(X)", "ssq(T)"),
                  ifelse(perc, "ssq(U)/ssq(Y)", "ssq(U)"), ifelse(perc,"sigH^2/ssq(U)","sigH^2"), "log LR",
                  "#steps", "last incr")
  print(round(outp,digits))
  return(invisible(outp))
}

#' Plot function for class PPLS.
#'
#' @param x A PPLS fit (an object of class PPLS)
#' @param XorY Plot loadings for X or Y?
#' @param i Positive integer. The i-th loading/score will be plotted on the first axis.
#' @param j Positive integer or NULL (default). The j-th loading/score will be plotted on the second axis, or if NULL, the i-th loading is plotted against its index.
#' @param use_ggplot2 Logical. Use ggplot2?
#' @param ... For consistency
#' @return If ggplot2 is active, returns the ggplot object. Else NULL.
#'
#' @details This function plots one or two loadings/scores.
#'
#' @export
plot.PPLS <- function (x, XorY = c("X", "Y"), i = 1, j = NULL, use_ggplot2=TRUE,...)
{
  if(ncol(x$W) < i || (!is.null(j) && ncol(x$W) < j)){stop("i and j cannot exceed #components!")}
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
    return(invisible(plt))
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
#' @details If you want just the X or Y scores, subset the result i.e. \code{scores.PPLS(fit,X,Y)[1:nrow(X),]} or \code{scores.PPLS(fit,X,Y)[-(1:nrow(X)),]}.
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


#' Performs one EM step (use \link{PPLS} for fitting a PPLS model)
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
#' @param Ipopu indicator of which populations present
#' @return A list with updated values for\itemize{
#'    \item{W }{Matrix}
#'    \item{C }{Matrix}
#'    \item{B }{Matrix}
#'    \item{sighat }{Vector containing updated sigX and sigY}
#'    \item{siglathat }{Vector containing sigH and sigT}
#'  }
#' @details This function passes its arguments to EMstepC (a C++ function, see \link{Rcpp}), which returns the expected sufficient statistics.
#' The maximization is done afterwards in this function. This may become a full C++ function later.
#'
#' @export
meta_EMstep <- function(X , Y , W. = W, C. = C, Ipopu, params)
{
  stopifnot(nrow(X) == length(Ipopu))
  Ipopu = as.factor(Ipopu)
  N = cumsum(table(Ipopu))
  W. = as.vector(W.)
  C. = as.vector(C.)

  ret1 = lapply(1:length(N), function(i){
    Ni = c(0,N)
    popui = (1+Ni[i]):Ni[i+1]
    X = X[popui,]
    Y = Y[popui,]
    e = with(params[[i]],{
      B_T. = B_T[1]
      sigX. = sigX[1]
      sigY. = sigY[1]
      sigH. = sigH[1]
      sigT. = sigT[1]

      g = sigT.^2*B_T.^2 + sigH.^2
      Kw = sigT.^2 - sigT.^4*B_T.^2/sigY.^2 + sigT.^4*B_T.^2*g/(sigY.^2*(g+sigY.^2))
      Kc = g - sigT.^4*B_T.^2/sigX.^2 + sigT.^6*B_T.^2/(sigX.^2*(sigT.^2+sigX.^2))
      Kwc = sigT.^2*B_T./(sigX.^2*sigY.^2) - Kc*sigT.^2*B_T./(sigX.^2*sigY.^2*(Kc+sigY.^2)) -
        sigT.^4*B_T./(sigX.^2*sigY.^2*(sigT.^2+sigX.^2)) +
        Kc*sigT.^4*B_T./(sigX.^2*sigY.^2*(Kc+sigY.^2)*(sigT.^2+sigX.^2))
      c1 = Kw / (sigX.^2*(Kw + sigX.^2));
      c3 = Kc / (sigY.^2*(Kc + sigY.^2));
      c2 = Kwc;
      meta_Estep(W.,C.,B_T.,X,Y,sigX.,sigY.,sigH.,sigT.,c1,c2,c3)
    })
    e2 = meta_Mstep(e)
    return(e2)
  })

  ret1$W. = orth(rowSums(sapply(ret1,function(e) c(sign(crossprod(ret1[[1]]$Cxt,e$Cxt)))*e$Cxt)))
  ret1$C. = orth(rowSums(sapply(ret1[-which(names(ret1)=="W.")],function(e) c(sign(crossprod(ret1[[1]]$Cxt,e$Cxt)))*e$Cyu)))

  return(ret1)
}

#' Performs PPLS fit in one direction. (use \link{PPLS} for fitting a PPLS model)
#'
#' @param X Numeric matrix.
#' @param Y Numeric matrix.
#' @param EMsteps strictly positive integer. Denotes the maximum number of EM steps to make.
#' @param atol Double, convergence criterium for the log-likelihood
#' @param initialGuess A string. Choose "o2m", "random" or "equal" depending on what type of initial guess you want.
#' @param customGuess A list with named components W,C,B,sigE,sigF,sigH,sigT.
#' @param critfunc Function measuring the increment value, i.e. f(L_{i+1} - L_{i}). Usually f == I or f == abs.
#' @param constraints List. Ideally being constructed with \link{fconstraint}.
#' @return A list with estimates for\itemize{
#'    \item{W}{Matrix}
#'    \item{C}{Matrix}
#'    \item{B}{Matrix}
#'    \item{sig}{ Vector containing updated sigX, sigY, sigH and sigT}
#'    \item{logvalue}{ Vector with last log-likelihood(s)}
#'    \item{Last_increment}{ Double, value of last increment.}
#'    \item{Number_steps}{ Integer, number of steps needed to stop.}
#'  }
#' @details This function estimates loadings and variances for one direction.
#'
#' @export
meta_PPLSi <- function(X,Y,Ipopu,EMsteps=1e2,atol=1e-4,initialGuess=c("equal","o2m","random","custom"),customGuess=NULL,critfunc=function(x){x},
                  constraints=fconstraint())
{
  stopifnot(nrow(X) == length(Ipopu))
  Ipopu = as.factor(Ipopu)
  stopifnot(length(constraints) == 7,all(names(constraints) == c("W","C","B","sigE","sigF","sigH","sigT")))
  p = ncol(X)
  q = ncol(Y)
  N = nrow(X)
  initialGuess = match.arg(initialGuess)
  if(!is.null(customGuess)){initialGuess = "custom"}
  if(initialGuess == "o2m"){
    sim=o2m(X,Y,1,0,0)
    Wnw=sim$W.;Cnw=sim$C.;Bnw=sim$B_T.[1];
    signw=sqrt(c((ssq(X)-ssq(sim$Tt))/N/p,(ssq(Y)-ssq(sim$U))/N/q));
    siglatnw=sqrt(c(ssq(sim$U)-ssq(sim$Tt*Bnw),ssq(sim$Tt))/N)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  } else if(initialGuess == "random"){
    Wnw=orth(runif(p));Cnw=orth(runif(q)); Bnw = rchisq(1,1); siglatnw =  rchisq(2,100)/100; signw = rchisq(2,10)/100
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }else if(initialGuess == "equal"){
    Wnw=orth(rep(1,p));Cnw=orth(rep(1,q)); Bnw = 1; siglatnw =  c(1,1); signw = c(1/p,1/q)
    #ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }else if(initialGuess == "custom"){
    stopifnot('list' %in% class(customGuess))
    Wnw = customGuess$W;Cnw=customGuess$C;Bnw=customGuess$B;siglatnw=with(customGuess,c(sigH,sigT));signw=with(customGuess,c(sigE,sigF))
  }
  Wnw = with(constraints,if(is.numeric(W)) W else Wnw )
  Cnw = with(constraints,if(is.numeric(C)) C else Cnw )
  #Bnw = with(constraints,if(is.numeric(B)) B else Bnw )
  #signw = with(constraints,c(if(is.numeric(sigE)) sigE else signw[1] , if(is.numeric(sigF)) sigF else signw[2]))
  #siglatnw = with(constraints,c(if(is.numeric(sigH)) sigH else siglatnw[1] , if(is.numeric(sigT)) sigT else siglatnw[2]))

  logvalue=matrix(NA, EMsteps+1, nlevels(Ipopu))

  logvalue[1,]=rep(logl_W(X,Y,Wnw,Cnw,Bnw,signw[1],signw[2],siglatnw[1],siglatnw[2]),nlevels(Ipopu))
  params = lapply(1:nlevels(Ipopu),function(i) list(B_T=Bnw,sigX=signw[1],sigY=signw[2],sigH=siglatnw[1],sigT=siglatnw[2]))
  Ni = c(0,cumsum(table(Ipopu)))
  datai = lapply(1:nlevels(Ipopu), function(j){
    popui = (1+Ni[j]):Ni[j+1]
    list(X = X[popui,], Y = Y[popui,])
  })
  for(i in 1:EMsteps){
#     if(any(signw < 100*.Machine$double.eps)){
#       return(list(W=NA,C=NA,B=NA,sig=NA,logvalue=NA,Last_increment = NA, Number_steps = i))
#     }
    fit=meta_EMstep(X,Y,Wnw,Cnw,Ipopu,params)
    params = lapply(1:nlevels(Ipopu), function(j){
      Bnw = (fit[[j]]$B)
      signw = fit[[j]]$sigh
      siglatnw = fit[[j]]$sigl

      #Wnw = with(constraints,if(is.numeric(W)) W else Wnw )
      #Cnw = with(constraints,if(is.numeric(C)) C else Cnw )
      #Bnw = with(constraints,if(is.numeric(B)) B else Bnw )
      #signw = with(constraints,c(if(is.numeric(sigE)) sigE else signw[1] , if(is.numeric(sigF)) sigF else signw[2]))
      #siglatnw = with(constraints,c(if(is.numeric(sigH)) sigH else siglatnw[1] , if(is.numeric(sigT)) sigT else siglatnw[2]))
      list(B_T=Bnw,sigX=signw[1],sigY=signw[2],sigH=siglatnw[1],sigT=siglatnw[2])
    })
    Wnw = (fit$W)
    #PYonw = (fit$PY)
    Cnw = (fit$C)
    #PXonw = (fit$PX)
    logvalue[i+1,] =
      sapply(1:nlevels(Ipopu),
             function(j) with(params[[j]],logl_W(datai[[j]]$X,datai[[j]]$Y,Wnw,Cnw,B_T,sigX,sigY,sigH,sigT)))
    if(critfunc(sum(logvalue[i+1,])-sum(logvalue[i,]))<atol){
      message("stopped at ",i,". Last incr ");print(logvalue[i+1,]-logvalue[i,]);
      break
      }
    #on.exit(return(list( W=c(Wnw),C=c(Cnw),B=Bnw,sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = logvalue[i]-logvalue[i-1], Number_steps = i)))
  }
  #vars = list(W = )

  #last_incr = logvalue[i+1]-logvalue[i]
  #if(any(diff(logvalue[0:i+1])<0)){warning("Not monotone")}
  #if(critfunc(logvalue[i+1]-logvalue[i])>=atol){warning(paste("Not converged, last increment was",critfunc(logvalue[i+1]-logvalue[i])))}
  #ret2 = list( W=c(Wnw),C=c(Cnw),B=c(Bnw),sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = last_incr, Number_steps = i)
  ret2 = list(W=c(Wnw),C=c(Cnw), params = params, log = logvalue[1:i+1,])
  return(ret2)
}

#' Performs the E step
#' Expectation step
#'
#' @param X First dataset.
#' @param Y Second dataset.
#' @param W X loadings, p times r.
#' @param C Y loadings, q times r.
#' @param B Diagonal regression matrix with positive elements.
#' @param sigE Positive number. Standard deviation of noise in X
#' @param sigF Positive number. Standard deviation of noise in Y
#' @param sigH Positive number. Standard deviation of noise in latent space in Y
#' @param sigT Diagonal positive matrix. Standard deviations(!) of latent space in X
#'
#' @return List with expected first and second moments of the latent variables and noise.
#' No aggregation is done yet, this means that all Cxx are matrices.
#' This may be useful to check assumptions of no correlation between the latent variables.
#' @export
Expect_M <- function(X,Y,W,C,B,sigE,sigF,sigH,sigT){
  N = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  a = ncol(W)

  covT = rbind(W%*%sigT^2, C%*%B%*%sigT^2)
  covU = rbind(W%*%B%*%sigT^2, C%*%B^2%*%sigT^2+sigH^2*C)
  invS= solve(sseXY_W(W,C,B,sigE,sigF,sigH,sigT))

  mu_T = cbind(X,Y) %*% invS %*% covT
  mu_U = cbind(X,Y) %*% invS %*% covU

  sigU = sqrt(sigT^2%*%B^2 + diag(sigH^2,ncol(C)))
  Ctt = sigT^2 - t(covT) %*% invS %*% covT + crossprod(mu_T) / N
  Cuu = sigU^2 - t(covU) %*% invS %*% covU + crossprod(mu_U) / N
  Cut = sigT^2 %*% B - t(covU) %*% invS %*% covT + crossprod(mu_U,mu_T) / N

  covE = rbind(diag(sigE^2,p), diag(0,q,p))
  mu_E = cbind(X,Y) %*% invS %*% covE
  Cee = diag(sigE^2,p) - t(covE) %*% invS %*% covE + crossprod(mu_E) / N

  covF = rbind(diag(0,p,q), diag(sigF^2,q))
  mu_F = cbind(X,Y) %*% invS %*% covF
  Cff = diag(sigF^2,q) - t(covF) %*% invS %*% covF + crossprod(mu_F) / N

  covH = rbind(0*W, sigH^2*C)
  mu_H = cbind(X,Y) %*% invS %*% covH
  Chh = diag(sigH^2,ncol(C)) - t(covH) %*% invS %*% covH + crossprod(mu_H) / N

  # mu_T = mu_U = Ctt = Cuu = Cut = Cee = Cff = Chh = NULL
  # for(i in 1:a){
  #   Efit = EMstep_W(X, Y, W[,i], C[,i], B[i,i], sigE, sigF, sigH, sigT[i,i])
  #   mu_T = cbind(mu_T,Efit$mu_T)
  #   mu_U = cbind(mu_U,Efit$mu_U)
  #   Ctt = c(Ctt,Efit$Ctt)
  #   Cut = c(Cut,Efit$Cut)
  #   Cuu = c(Cuu,Efit$Cuu)
  #   Cee = as.matrix(Efit$Cee)
  #   Cff = as.matrix(Efit$Cff)
  #   Chh = as.matrix(Efit$Chh)
  # }
  # Ctt = diag(Ctt,a)
  # Cut = diag(Cut,a)
  # Cuu = diag(Cuu,a)

  list(mu_T = mu_T, mu_U = mu_U, Ctt = Ctt, Cuu = Cuu,
       Cut = Cut, Cee = Cee, Cff = Cff, Chh = Chh)
}

#' The M step
#'
#' Maximization step
#'
#' @inheritParams Expect_M
#' @param fit A list as produced by \code{\link{Expect_M}}
#'
#' @return A list with updated estimates W, C, B, sigE, sigF, sigH and sigT.
#'
#' @export
Maximiz_M <- function(fit,X,Y, type = "SVD"){
  outp = with(fit,{
    list(
      W = orth(t(X) %*% mu_T,type=type),
      C = orth(t(Y) %*% mu_U,type=type),
      B = Cut %*% solve(Ctt) * diag(1,nrow(Cut)),
      sigE = sqrt(tr(Cee)/ncol(Cee)),
      sigF = sqrt(tr(Cff)/ncol(Cff)),
      sigH = sqrt(tr(Chh)/ncol(Chh)),
      sigT = sqrt(Ctt * diag(1,nrow(Ctt)))
    )
  })
  return(outp)
}

#' The PPLS fitting function
#'
#' Simultaneous PPLS fitting function
#'
#' @inheritParams Expect_M
#' @inheritParams Maximiz_M
#' @param a Positive integer, number of components
#' @param EMsteps Positive integer, number of EM steps
#' @param atol positive double, convergence criterium
#'
#' @return list of class PPLS_simul with expectations, loglikelihoods and estimates.
#'
#' @export
PPLS_simult <- function(X, Y, a, EMsteps = 10, atol = 1e-4, type = "SVD"){
  p = ncol(X)
  q = ncol(Y)
  W. = orth(matrix(rnorm(p*a),p)) #svd(X,nu=0,nv=a)$v
  C. = orth(matrix(rnorm(q*a),q)) #svd(Y,nu=0,nv=a)$v
  B. = diag(sort(rnorm(a,mean = 1, sd = .5)),a)
  sigE. = 1/p
  sigF. = 1/q
  sigH. = 0.1/a
  sigT. = diag(1,a)
  logl_incr = 1:EMsteps*NA
  for(i in 1:EMsteps){
    Expect_M(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.) %>% Maximiz_M(X,Y) -> outp
    W. = outp$W
    C. = outp$C
    B. = outp$B
    sigE. = outp$sigE
    sigF. = outp$sigF
    sigH. = outp$sigH
    sigT. = outp$sigT

    logl_incr[i] = logl_W(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.)
    if(i > 1 && diff(logl_incr)[i-1] < atol){ break}
  }
  signLoad = sign(diag(sigT. %*% B.))
  rotLoad = order(diag(sigT. %*% B. %*% diag(signLoad)), decreasing=TRUE)
  outp$W = W.[,rotLoad] %*% diag(signLoad, a)
  outp$C = C.[,rotLoad] %*% diag(signLoad, a)
  outp$B = diag(diag(B. %*% diag(signLoad, a))[rotLoad])
  outp$sigT = diag(diag(sigT.)[rotLoad])
  logl_incr = logl_incr[1:i]
  if(any(diff(logl_incr) < 0)) warning("Negative increments of likelihood")
  Eout = Expect_M(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.)
  outpt = list(Expectations = Eout, loglik = logl_incr, estimates = outp)

  class(outpt) <- "PPLS_simult"
  outpt
}

# mseLoadings <- function(W0, W1, W2 = NULL, W3 = NULL){
#   a = ncol(W0)
#   W1 = W1%*%sign(crossprod(W1,W0)*diag(a))
#   if(!is.null(W2)) W2 = W2%*%sign(crossprod(W2,W0)*diag(a))
#   if(!is.null(W3)) W3 = W3%*%sign(crossprod(W3,W0)*diag(a))
#
#   sapply(list(W1 = W1, W2 = W2, W3 = W3),function(e){if(!is.null(e)) mse(W0, e)})
# }

#' Variances for PPLS_simult
#'
#' Calculates asymptotic variances for PPLS loadings
#'
#' @inheritParams Expect_M
#' @inheritParams Maximiz_M
#' @inheritParams PPLS_simult
#' @param data Data matrix, X or Y
#' @param XorY Which data matrix did you supply, X or Y?
#'
#' @return SE's for the loadings of corresponding data (so W or C)
variances.PPLS_simult <- function(fit, data, XorY = c("X", "Y")){
  W = fit$estim[ifelse(XorY=="X", "W", "C")][[1]]
  X = data
  Ctt = fit$Expec[ifelse(XorY=="X", "Ctt", "Cuu")][[1]]
  mu_T = fit$Expec[ifelse(XorY=="X", "mu_T", "mu_U")][[1]]
  outp <- lapply(1:a, function(i) {
    W = as.matrix(W[,i])
    Ctt = t(Ctt[i,i])
    mu_T = mu_T[,i]
    Vt = Ctt - crossprod(mu_T)
    Cxt = t(X) %*% mu_T
    B_star = c(Ctt)/(4*fit$estim$sigE^2) * (diag(1,ncol(X)) - tcrossprod(W))

    SSt_expec = t(X) %*% diag(c(Ctt),nrow(X)) %*% X - t(X) %*% mu_T %*% (Ctt + 2*Vt) %*% t(W) -
      W %*% (Ctt + 2*Vt) %*% t(mu_T) %*% X + W %*% (Ctt^2 + 4*t(mu_T)%*%mu_T +2*Vt) %*% t(W)
    SSt_expec = SSt_expec/(4*fit$estim$sigE^4)

    SSt_star = tcrossprod(Cxt - W*c(Ctt))#Cxt%*%t(Cxt) - Cxt%*%Ctt%*%t(W) - W%*%Ctt%*%t(Cxt) + W%*%Ctt^2%*%t(W)
    SSt_star = SSt_star/(4*fit$estim$sigE^4)

    Iobs = B_star - SSt_expec + SSt_star
    list(B_exp = B_star, SSt_exp = SSt_expec, SSt_star = SSt_star)
  })
  outp$seLoad = sapply(1:a, function(i) with(outp[[i]], sqrt(-diag(solve(B_exp - SSt_exp + SSt_star)))))
  outp
  #aparte matrices
}

