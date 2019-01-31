source('https://raw.githubusercontent.com/selbouhaddani/PO2PLS/master/R/PO2PLS_functions.R')
tr <- function(X) sum(diag(X))
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
Expect_M_DA <- function(X,Y,W,C,B,sigE,sigF,sigH,sigT,debug=T){
  N = nrow(X)
  p = ncol(X)
  q = ncol(as.matrix(Y))
  a = ncol(W)

  if(debug){
    covT = rbind(W%*%sigT^2, C%*%t(B)%*%sigT^2)
    covU = rbind(W%*%sigT^2%*%B, C%*%t(B)%*%sigT^2%*%B+sigH^2*C)
    SX. = tcrossprod(W%*%sigT)+sigE^2*diag(1,nrow=p)
    SXY. = W%*%sigT^2%*%B%*%t(C)
    SY. = tcrossprod(C%*%t(B)%*%sigT)+tcrossprod(C)*sigH^2+sigF^2*diag(1,nrow=q)
    #invSX. = woodinv(sigX.,cbind(W.,P_Yosc.))
    #invSY. = woodinv(sigY.,cbind(C.,P_Xosc.))
    Sfull = blockm(SX.,SXY.,SY.)
    invS= solve(Sfull)

    mu_T = cbind(X,Y) %*% invS %*% covT
    mu_U = cbind(X,Y) %*% invS %*% covU

    sigU = sqrt(abs(t(B)%*%sigT^2%*%B + sigH^2))
    Ctt = sigT^2 - t(covT) %*% invS %*% covT + crossprod(mu_T) / N
    Cuu = sigU^2 - t(covU) %*% invS %*% covU + crossprod(mu_U) / N
    Cut = t(B) %*% sigT^2 - t(covU) %*% invS %*% covT + crossprod(mu_U,mu_T) / N

    covE = rbind(diag(sigE^2,p), diag(0,q,p))
    mu_E = cbind(X,Y) %*% invS %*% covE
    Cee = tr(diag(sigE^2,p) - t(covE) %*% invS %*% covE + crossprod(mu_E) / N)/p

    covF = rbind(diag(0,p,q), diag(sigF^2,q))
    mu_F = cbind(X,Y) %*% invS %*% covF
    Cff = tr(diag(sigF^2,q) - t(covF) %*% invS %*% covF + crossprod(mu_F) / N)/q

    covH = rbind(matrix(0,p,q), sigH^2*C)
    mu_H = cbind(X,Y) %*% invS %*% covH
    Chh = sigH^2 - t(covH) %*% invS %*% covH + crossprod(mu_H) / N

  } else {
    ##############
    sigT = diag(sigT)
    B = diag(B)
    g = sapply(1:a, function(i) sigT[i]^2*B[i]^2 + sigH^2)
    Kw = sapply(1:a, function(i) sigT[i]^2 - sigT[i]^4*B[i]^2/sigF^2 + sigT[i]^4*B[i]^2*g[i]/(sigF^2*(g[i]+sigF^2)))
    Kc = sapply(1:a, function(i) g[i] - sigT[i]^4*B[i]^2/sigE^2 + sigT[i]^6*B[i]^2/(sigE^2*(sigT[i]^2+sigE^2)))
    Kwc = sapply(1:a, function(i) sigT[i]^2*B[i]/(sigE^2*sigF^2) - Kc[i]*sigT[i]^2*B[i]/(sigE^2*sigF^2*(Kc[i]+sigF^2)) -
                   sigT[i]^4*B[i]/(sigE^2*sigF^2*(sigT[i]^2+sigE^2)) +
                   Kc[i]*sigT[i]^4*B[i]/(sigE^2*sigF^2*(Kc[i]+sigF^2)*(sigT[i]^2+sigE^2)))
    c1 = sapply(1:a, function(i) Kw[i] / (sigE^2*(Kw[i] + sigE^2)))
    c3 = sapply(1:a, function(i) Kc[i] / (sigF^2*(Kc[i] + sigF^2)))
    c2 = Kwc
    #    if(global_debug) return(c(c1, c2, c3))
    sigT = diag(sigT, a)
    B = diag(B, a)
    c1 = diag(c1, a)
    c2 = diag(c2, a)
    c3 = diag(c3, a)

    varU = sigT^2%*%B^2 + diag(sigH^2,a)
    Xw = X%*%W
    Yc = Y%*%C
    mu_T = sigE^-2 * Xw%*%sigT^2 + sigF^-2 * Yc%*%sigT^2%*%B - Xw%*%c1%*%sigT^2 -
      Xw%*%c2%*%sigT^2%*%B - Yc%*%c2%*%sigT^2 - Yc%*%c3%*%B%*%sigT^2
    mu_U = sigE^-2 * Xw%*%sigT^2%*%B + sigF^-2 * Yc%*%varU -
      Xw%*%c1%*%sigT^2%*%B - Xw%*%c2%*%varU - Yc%*%c2%*%sigT^2%*%B - Yc%*%c3%*%varU

    Ctt = sigT^2 - sigE^-2*sigT^4 - sigF^-2*sigT^4%*%B^2 + sigT^4%*%c1 + 2*sigT^4%*%B%*%c2 +
      sigT^4%*%B^2%*%c3 + crossprod(mu_T) / N
    Cuu = varU - sigE^-2*sigT^4%*%B^2 - sigF^-2*varU^2 + sigT^4%*%B^2%*%c1 +
      2*sigT^2%*%B%*%varU%*%c2 + varU^2%*%c3 + crossprod(mu_U) / N
    Cut = sigT^2 %*% B - sigE^-2*sigT^4%*%B - sigF^-2*sigT^2%*%B%*%varU + sigT^4%*%B%*%c1 +
      sigT^2%*%varU%*%c2 + sigT^4%*%B^2%*%c2 + sigT^2%*%B%*%varU%*%c3 + crossprod(mu_U,mu_T) / N

    mu_E = X - sigE^2*Xw%*%c1%*%t(W) - sigE^2*Yc%*%c2%*%t(W)
    #   ssq_E = ssq(X) - 2*sigE^2*(tr(crossprod(Xw, Xw%*%c1)) - tr(crossprod(Xw, Yc%*%c2))) +
    #     sigE^4*(2*tr(crossprod(Xw, Yc%*%c1%*%c2)) + ssq(Xw%*%c1) + ssq(Yc%*%c2))
    Cee = (p*sigE^2 - p*sigE^2 + sigE^4*sum(c1) + ssq(mu_E) / N) / p

    mu_F = Y - sigF^2*Yc%*%c3%*%t(C) - sigF^2*Xw%*%c2%*%t(C)
    Cff = (q*sigF^2 - q*sigF^2 + sigF^4*sum(c3) + ssq(mu_F) / N) / q

    mu_H = sigF^-2*sigH^2*Yc - sigH^2*(Xw%*%c2 + Yc%*%c3)
    Chh = diag(sigH^2 - sigH^4/sigF^2,a) + sigH^4*c3 + crossprod(mu_H) / N
    ##############
  }
  list(mu_T = mu_T, mu_U = mu_U, Ctt = (Ctt), Cuu = (Cuu),
       Cut = Cut, Cee = as.matrix(Cee), Cff = as.matrix(Cff), Chh = (Chh))
}

#' The M step
#'
#' Maximization step
#'
#' @inheritParams Expect_M
#' @param fit A list as produced by \code{\link{Expect_M}}
#' @param type String. One of "SVD" or "QR"
#' @return A list with updated estimates W, C, B, sigE, sigF, sigH and sigT.
#'
#' @export
Maximiz_M_DA <- function(fit,X,Y, type = c("SVD","QR")){
  type <- match.arg(type)
  outp = with(fit,{
    list(
      W = orth(t(X) %*% mu_T,type=type),
      C = orth(t(Y) %*% mu_U,type=type),
      B = Cut %*% solve(Ctt) %>% t,
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
#' @param ... debug, logical. Should the slow Expect_M be used?
#'
#' @return list of class PPLS_simul with expectations, loglikelihoods and estimates.
#'
#' @export
PPLS_simult_DA <- function(X, Y, a, EMsteps = 10, atol = 1e-4, type = c("SVD","QR"), ...){
  p = ncol(X)
  Y = as.matrix(Y)
  q = ncol(as.matrix(Y))
  type = match.arg(type)
  init_fit <- svd(X,nu=a,nv=a)
  W. = init_fit$v # orth(matrix(rnorm(p*a),p)) #svd(X,nu=0,nv=a)$v
  C. = diag(1,ncol(Y))  # orth(matrix(rnorm(q*a),q)) #svd(Y,nu=0,nv=a)$v
  B. = solve(cov(init_fit$u))%*%cov(init_fit$u,Y) # diag(sort(abs(rnorm(a,0, .5)),T),a)
  sigE. = sum(init_fit$d[-(1:a)]^2)%>%sqrt
  sigF. = 0 #f0$sig[a,2] # 1/q
  sigH. = sd(Y - init_fit$u%*%B.)
  sigT. = diag(init_fit$d[1:a],a) # diag(1,a)

#  logl_incr = 1:EMsteps*NA
  err <- 1:EMsteps*NA
  for(i in 1:EMsteps){
    Expect_M_DA(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.,...) %>% Maximiz_M_DA(X,Y,type) -> outp
    err[i] <- crossprod(outp$W, W.) %>% abs %>% max
    W. = outp$W
    C. = outp$C
    B. = outp$B
    sigE. = outp$sigE
    sigF. = outp$sigF
    sigH. = outp$sigH
    sigT. = outp$sigT
    #logl_incr[i] = logl_W(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.)
    #if(i > 1 && diff(logl_incr)[i-1] < atol){ break}
    #if(i > 1 && diff(err)[i-1] < (1 - atol) ){ break }
  }
  # signLoad = sign(diag(sigT. %*% B.))
  # rotLoad = order(diag(sigT. %*% B. %*% diag(signLoad,a)), decreasing=TRUE)
  # outp$W = W.[,rotLoad] %*% diag(signLoad, a)
  # outp$C = C.[,rotLoad] %*% diag(signLoad, a)
  # outp$B = diag(diag(B. %*% diag(signLoad, a))[rotLoad],a)
  # outp$sigT = diag(diag(sigT.)[rotLoad],a)
  # logl_incr = logl_incr[1:i]
  # if(any(diff(logl_incr) < 0)) warning("Negative increments of likelihood")
  # Eout = Expect_M(X,Y,W.,C.,B.,sigE.,sigF.,sigH.,sigT.)
  outpt = list(estimates = outp, error = err[1:i])
  #
  # class(outpt) <- "PPLS_simult"
  outpt
}

library(magrittr)
library(OmicsPLS)
library(PPLS)
N <- 50
N_test <- 1000
p <- 20
q <- 3
r <- 2
Tt <- matrix(rnorm(N*r),N)
W = orth(matrix(rt(p*r,2),p))
X <- tcrossprod(Tt, W)
X <- X + matrix(rnorm(prod(dim(X)), sd=0.2*sd(X)),nrow(X))
X %<>% scale(scale=T)
B <- diag(1,r,q)
Y <- Tt %*% B
Y <- Y + matrix(rnorm(prod(dim(Y)),sd = 0.2*sd(Y)),nrow(Y))
Y %<>% scale(scale=T)

#X <- matrix(rnorm(N*p),N) %>% scale
#Y <- X %*% orth(cbind(1:p, rep(1,p))) %>% scale

fit <- PPLS_simult_DA(X,Y,r,1e3,1e-6)
fit2 <- PPLS_simult(X, Y, r, 1e3, 1e-6)
#crossprod(fit$est$W, W)


Tt_test <- matrix(rnorm(N_test*r),N_test)
X_test <- tcrossprod(Tt_test, W)
X_test <- X_test + matrix(rnorm(prod(dim(X_test)), sd=0.2*sd(X_test)),nrow(X_test))
X_test %<>% scale(scale=T)
Y_test <- Tt_test %*% B
Y_test <- Y_test + matrix(rnorm(prod(dim(Y_test)),sd = 0.2*sd(Y_test)),nrow(Y_test))
Y_test %<>% scale(scale=T)

#X_test <- matrix(rnorm(N_test*p),N_test) %>% scale
#Y_test <- X_test %*% orth(cbind(1:p, rep(1,p))) %>% scale

ssq(Y_test - predict(pls::plsr(Y~X, ncomp=r), X_test, ncomp=r)[,,1])/ssq(Y_test)
ssq(Y_test - X_test %*% coef(lm(Y~X-1)))/ssq(Y_test)
ssq(Y_test - X_test %*% fit$es$W %*% fit$es$B)/ssq(Y_test)
ssq(Y_test - X_test %*% fit2$es$W %*% fit2$es$B %*% t(fit2$es$C))/ssq(Y_test)

