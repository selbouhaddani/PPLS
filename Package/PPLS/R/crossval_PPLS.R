#' Predicts X or Y with PPLS
#'
#' Predicts X or Y based on new data on Y or X
#'
#' @param object PPLS fit
#' @param newdata New data, which one of X or Y is specified in \code{XorY}.
#' @param XorY Character specifying whether \code{newdata} is X or Y.
#'
#' @return Predicted Data
#'
#' @export
predict.PPLS <- function(object, newdata, XorY = c("X","Y"), ...) {
  XorY = match.arg(XorY)
  if(!is.matrix(newdata)) newdata <- t(newdata)
  switch(XorY,
         X = if(ncol(newdata) != nrow(object$W)) stop("Number of columns mismatch!"),
         Y = if(ncol(newdata) != nrow(object$C)) stop("Number of columns mismatch!"))

  pred = switch(XorY,
                Y = with(object, (newdata) %*% C %*% diag(1/B,length(B)) %*% t(W)),
                X = with(object, (newdata) %*% W %*% diag(B,length(B)) %*% t(C)))

  return(pred)
}

#' Predicts X or Y with PPLS
#'
#' Predicts X or Y based on new data on Y or X
#'
#' @inheritParams predict.PPLS
#' @inheritParams PPLS
#' @param X X data
#' @param Y Y data
#' @param nr_folds Positive integer. Number of folds to consider. Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two folds does not make sense.
#' @param ... PPLS arguments
#'
#' @return Cross-validated RMSEP averaged over all elements
#'
#' @export
cv_PPLS <- function(X, Y, nr_comp, nr_folds, ...){
  kcv <- nr_folds
  N <- nrow(X)
  blocks <- as.numeric(cut(seq(1:N), breaks=kcv, labels=F))
  folds <- sample(N)
  err <- NA*1:kcv
  for(i in 1:kcv){
    ii <- which(blocks==i)
    fit <- PPLS(X[-folds[ii],], Y[-folds[ii],], nr_comp, ...)
    err[i] <- mean(c(sqrt(mse(Y[folds[ii],], predict(fit, X[folds[ii],], "X")))))
  }
  mean(err)
}


#' Cross-validate procedure for PPLS
#'
#' @inheritParams PPLS
#' @inheritParams cv_PPLS
#' @param a Vector of positive integers. Denotes the numbers of joint components to consider.
#' @param nr_cores Positive integer. Number of cores to use for CV. You might want to use \code{\link{detectCores}()}. Defaults to 1.
#'
#' @details This is the standard CV approach. It minimizes the sum of the prediction errors of X and Y over a three-dimensional grid of integers.
#' Parallelization is possible on all platforms. On Windows it uses \code{\link{makePSOCKcluster}}, then exports all necessary objects to the workers,
#'  and then calls \code{\link{parLapply}}. On OSX and Linux the more friendly \code{\link{mclapply}} is used, which uses forking.
#'  A print method is defined, printing the minimizers and minimum in a readible way. Also the elapsed time is tracked and reported.
#'
#' @return List of class \code{"cvo2m"} with the original and sorted Prediction errors and the number of folds used.
#'
#' @examples
#' local({
#' X = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' Y = scale(jitter(tcrossprod(rnorm(100),runif(10))))
#' crossval_PPLS(X, Y, a = 1:4,
#'              nr_folds = 5, nr_cores = 1, EMsteps = 1e3)
#' })
#'
#' @export
crossval_PPLS <- function(X, Y, a, nr_folds, nr_cores,...) {
  tic = proc.time()
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...")}
  kcv = nr_folds
  ax = 0
  ay = 0
  stopifnot(ncol(X) > max(a)+max(ax) , ncol(Y) > max(a)+max(ay) , nrow(X) >= kcv)
  stopifnot(nr_cores == abs(round(nr_cores)))
  if(nr_folds==1){stop("Cross-validation with 1 fold does not make sense, use 2 folds or more")}

  #parms = data.frame(nx = ax)
  #parms = merge(parms,data.frame(ny = ay))
  parms = data.frame(a = a)
  parms = apply(parms,1,as.list)
  cl_crossval_PPLS <- NULL

  on.exit({if(!is.null(cl_crossval_PPLS)) stopCluster(cl_crossval_PPLS)})

  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_PPLS <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_PPLS, library(PPLS))
    clusterExport(cl_crossval_PPLS, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_PPLS,parms,function(e){
      cv_PPLS(X, Y, e$a, nr_folds, ...)
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      cv_PPLS(X, Y, e$a, nr_folds, ...)
    })
  }
  outp = list(errors = unlist(outp), which_error = which.min(outp))
  class(outp) <- "cvPPLS"
  toc = proc.time() - tic
  outp$time = round(toc[3],2)
  message("Min at a = ", outp$which_error)
  return(outp)
}
