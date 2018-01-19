summary.o2m <- function(object, digits = 3, ...) {
  fit <- object
  a <- ncol(fit$W.)
  if(digits != round(digits) || digits <= 0) stop("digits must be a positive integer")
  outp <- with( fit, list(
    Comp = a,
    R2_X = R2X,
    R2_Y = R2Y,
    R2_Xjoint = R2Xcorr,
    R2_Yjoint = R2Ycorr,
    R2_Xhat = R2Xhat,
    R2_Yhat = R2Yhat,
    R2_Xpred = R2Xhat / R2Xcorr,
    R2_Ypred = R2Yhat / R2Ycorr,
    B_T = B_T.,
    B_U = B_U,
    flags = flags,
    digits = digits
  ) )
  class(outp) <- "summary.o2m"
  #   Mname <- list(c(""), c("Comp", "R2X", "R2Y", "R2Xcorr", "R2Ycorr", "R2Xhat", "R2Yhat", "XRatio", "YRatio"))
  #   M <- matrix(c(ifelse(perc,a/100,a), fit$R2X, fit$R2Y, fit$R2Xcorr, fit$R2Ycorr, fit$R2Xhat, fit$R2Yhat, fit$R2Xhat/fit$R2Xcorr, 
  #                 fit$R2Yhat/fit$R2Ycorr), nrow = 1, dimnames = Mname)
  #   return(round((1 + perc * 99) * M, 4 - perc * 2))
  outp
}

PPLS_to_o2m <- function(X_true, Y_true, fit_PPLS, stripped = TRUE){
  
  Xnames = dimnames(X_true)
  Ynames = dimnames(Y_true)
  
  W <- fit_PPLS$W
  C <- fit_PPLS$C
  B_T <- diag(fit_PPLS$B, length(fit_PPLS$B))
  B_U <- solve(B_T)
  Tt <- X_true %*% W
  U <- Y_true %*% C
  T_Yosc <- U_Xosc <- matrix(0, nrow(Tt), 1)
  P_Yosc <- W_Yosc <- matrix(0, nrow(W), 1)
  P_Xosc <- C_Xosc <- matrix(0, nrow(C), 1)
  
  H_TU <- 0*Tt# - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  ssqX <- ssq(X_true)
  ssqY <- ssq(Y_true)
  
  R2Xcorr <- ssq(Tt) / ssqX
  R2Ycorr <- ssq(U) / ssqY
  R2X_YO <- ssq(T_Yosc %*% t(P_Yosc)) / ssqX
  R2Y_XO <- ssq(U_Xosc %*% t(P_Xosc)) / ssqY
  R2Xhat <- (ssq(U %*% B_U %*% t(W)) / ssqX)
  R2Yhat <- (ssq(Tt %*% B_T %*% t(C)) / ssqY)
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
                T_Yosc. = T_Yosc, U_Xosc. = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
                B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
                R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- c("o2m","o2m_stripped")
  
  #toc <- proc.time() - tic
  model$flags = c(time = NA, 
                  list(n = ncol(W), nx = 0, ny = 0, 
                       stripped = TRUE, highd = FALSE, 
                       call = match.call(), ssqX = ssqX, ssqY = ssqY,
                       varXjoint = apply(model$Tt,2,ssq),
                       varYjoint = apply(model$U,2,ssq),
                       varXorth = apply(model$P_Y,2,ssq)*apply(model$T_Y,2,ssq),
                       varYorth = apply(model$P_X,2,ssq)*apply(model$U_X,2,ssq)))
  return(model)
}

PPLS_simult_to_o2m <- function(X_true, Y_true, fit_PPLS, stripped = TRUE){
  
  Xnames = dimnames(X_true)
  Ynames = dimnames(Y_true)
  
  N <- nrow(X_true)
  p <- ncol(X_true)
  q <- ncol(Y_true)
  r <- ncol(fit_PPLS$est$W)
  
  trc <- function(X){ sum(diag(as.matrix(X))) }
  
  W <- fit_PPLS$est$W
  C <- fit_PPLS$est$C
  B_T <- fit_PPLS$est$B
  B_U <- solve(B_T)
  Tt <- fit_PPLS$Exp$mu_T
  U <- fit_PPLS$Exp$mu_U
  T_Yosc <- U_Xosc <- matrix(0, nrow(Tt), 1)
  P_Yosc <- W_Yosc <- matrix(0, nrow(W), 1)
  P_Xosc <- C_Xosc <- matrix(0, nrow(C), 1)
  
  H_TU <- 0*Tt# - U %*% B_U
  H_UT <- U - Tt %*% B_T
  
  ssqX <- ssq(X_true)
  ssqY <- ssq(Y_true)
  
  R2Xcorr <- with(fit_PPLS$est, ssq(sigT^2) / (ssq(sigT^2) + p*sigE^2))
  R2Ycorr <- with(fit_PPLS$est, ssq(sigT^2 %*% B^2 + diag(sigH^2,r)) / (ssq(sigT^2 %*% B^2 + diag(sigH^2,r)) + q*sigF^2))
  R2X_YO <- 0
  R2Y_XO <- 0
  R2Xhat <- NA
  R2Yhat <- with(fit_PPLS$est, ssq(sigT^2 %*% B) / (ssq(sigT^2 %*% B^2) + r*sigH^2 + q*sigF^2))
  R2X <- R2Xcorr + R2X_YO
  R2Y <- R2Ycorr + R2Y_XO
  
  rownames(Tt) <- rownames(T_Yosc) <- rownames(H_TU) <- Xnames[[1]]
  rownames(U) <- rownames(U_Xosc) <- rownames(H_TU) <- Ynames[[1]]
  rownames(W) <- rownames(P_Yosc) <- Xnames[[2]]
  rownames(C) <- rownames(P_Xosc) <- Ynames[[2]]
  
  model <- list(Tt = Tt, U = U, W. = W, C. = C, P_Yosc. = P_Yosc, P_Xosc. = P_Xosc,
                T_Yosc. = T_Yosc, U_Xosc. = U_Xosc, W_Yosc = W_Yosc, C_Xosc = C_Xosc,
                B_T. = B_T, B_U = B_U, H_TU = H_TU, H_UT = H_UT, 
                R2X = R2X, R2Y = R2Y, R2Xcorr = R2Xcorr, R2Ycorr = R2Ycorr, 
                R2Xhat = R2Xhat, R2Yhat = R2Yhat)
  class(model) <- c("o2m","o2m_stripped")
  
  #toc <- proc.time() - tic
  model$flags = c(time = NA, 
                  list(n = ncol(W), nx = 0, ny = 0, 
                       stripped = TRUE, highd = FALSE, 
                       call = match.call(), ssqX = ssqX, ssqY = ssqY,
                       varXjoint = apply(model$Tt,2,ssq),
                       varYjoint = apply(model$U,2,ssq),
                       varXorth = apply(model$P_Y,2,ssq)*apply(model$T_Y,2,ssq),
                       varYorth = apply(model$P_X,2,ssq)*apply(model$U_X,2,ssq)))
  return(model)
}