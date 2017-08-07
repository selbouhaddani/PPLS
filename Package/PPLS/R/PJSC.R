tr <- function(X){
  if(!is.matrix(X)) X = as.matrix(X)
  if(nrow(X) != ncol(X)) warning("X not symmetric")
  sum(diag(X))
}

blockm<-function(A,B,C)
  #input: Matrices A,B,C
  #output: the block matrix
  # A    B
  #t(B)  C
{
  M = rbind(cbind(A,B),cbind(t(B),C))
  return(M)
}

l_step = function(X, Y, W, C, Wo, Co, Phi1, Phi2, Psi1, Psi2) {
  N = nrow(X)
  W = as.matrix(W)
  Wo = as.matrix(Wo)
  C = as.matrix(C)
  Co = as.matrix(Co)
  r = ncol(W)
  p = ncol(X)
  q = ncol(Y)

  Sigma = blockm(
    W %*% (diag(1, r) + Phi1) %*% t(W) + Wo %*% t(Wo) + Psi1,
    W %*% t(C),
    C %*% (diag(1, r) + Phi2) %*% t(C) + Co %*% t(Co) + Psi2
  )
  #S = crossprod(X,Y)/N
  invSigma = solve(Sigma)

  - N * (p + q) / 2 * log(2 * pi) - N / 2 * determinant(Sigma)[[1]] - 1 /
    2 * tr(crossprod(cbind(X, Y)) %*% invSigma)

}

cov_step = function(X, Y, W, C, Wo, Co, Phi1, Phi2, Psi1, Psi2){
  r = ncol(W)
  Sj = blockm(
    W %*% (diag(1, r) + Phi1) %*% t(W),
    W %*% t(C),
    C %*% (diag(1, r) + Phi2) %*% t(C)
  )
  So = blockm(
    tcrossprod(Wo),
    matrix(0,p,q),
    tcrossprod(Co)
  )
  Sr = blockm(
    Psi1,
    matrix(0,p,q),
    Psi2
  )
  list(joint=Sj, systematic=So, residual=Sr)
}

E_step <- function(X, Y, W, C, Wo, Co, Phi1, Phi2, Psi1, Psi2) {
  N = nrow(X)
  W = as.matrix(W)
  Wo = as.matrix(Wo)
  C = as.matrix(C)
  Co = as.matrix(Co)
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  p = ncol(X)
  q = ncol(Y)
  XY = cbind(X, Y)
  Sigma = blockm(
    W %*% (diag(1, r) + Phi1) %*% t(W) + Wo %*% t(Wo) + Psi1,
    W %*% t(C),
    C %*% (diag(1, r) + Phi2) %*% t(C) + Co %*% t(Co) + Psi2
  )
  #S = crossprod(X,Y)/N
  invSigma = solve(Sigma)
  #invS1 = solve(W %*% (diag(1,r) + Phi1) %*% t(W) + Wo %*% t(Wo) + Psi1)
  #invS2 = solve(C %*% (diag(1,r) + Phi2) %*% t(C) + Co %*% t(Co) + Psi2)

  ############ X part
  CovXT = rbind(W %*% (diag(1, r) + Phi1), C)
  muT = (XY) %*% invSigma %*% CovXT
  muTo = XY %*% invSigma %*% rbind(Wo, matrix(0, q, rx))

  Cxt = t(X) %*% muT / N
  Cxto = t(X) %*% muTo / N

  Ctt = (diag(1, r) + Phi1) -
    t(CovXT) %*% invSigma %*% CovXT + crossprod(muT) / N
  Ctoto = diag(1, rx) -
    t(rbind(Wo, matrix(0, q, rx))) %*% invSigma %*% rbind(Wo, matrix(0, q, rx)) +
    crossprod(muTo) / N

  ########### Y part
  CovYU = rbind(W, C %*% (diag(1, r) + Phi2))
  muU = (XY) %*% invSigma %*% CovYU
  muUo = XY %*% invSigma %*% rbind(matrix(0, p, ry), Co)

  Cyu = t(Y) %*% muU / N
  Cyuo = t(Y) %*% muUo / N

  Cuu = (diag(1, r) + Phi2) -
    t(CovYU) %*% invSigma %*% CovYU + crossprod(muU) / N
  Cuouo = diag(1, ry) -
    t(rbind(matrix(0, p, ry), Co)) %*% invSigma %*% rbind(matrix(0, p, ry), Co) +
    crossprod(muUo) / N

  ########## Noise part
  Psi12 = rbind(Psi1, matrix(0, q, p))
  Psi21 = rbind(matrix(0, p, q), Psi2)
  Phi12 = rbind(W %*% Phi1, matrix(0, q, r))
  Phi21 = rbind(matrix(0, p, r), C %*% Phi2)
  Cee = Psi1 - t(Psi12) %*% invSigma %*% Psi12 + t(Psi12) %*% invSigma %*% crossprod(XY) %*% invSigma %*% Psi12 / N
  Cff = Psi2 - t(Psi21) %*% invSigma %*% Psi21 + t(Psi21) %*% invSigma %*% crossprod(XY) %*% invSigma %*% Psi21 / N
  Cetet = Phi1 - t(Phi12) %*% invSigma %*% Phi12 + t(Phi12) %*% invSigma %*% crossprod(XY) %*% invSigma %*% Phi12 / N
  Ceueu = Phi2 - t(Phi21) %*% invSigma %*% Phi21 + t(Phi21) %*% invSigma %*% crossprod(XY) %*% invSigma %*% Phi21 / N

  #if(any(diag(Cetet)<0)) {warning("1");Cetet = - t(W) %*% invS1 %*% W + 2 * t(W) %*% invS1 %*% crossprod(X) %*% invS1 %*% W / N}
  #if(any(diag(Ceueu)<0)) {warning("2");Ceueu = - t(C) %*% invS2 %*% C + 2 * t(C) %*% invS2 %*% crossprod(Y) %*% invS2 %*% C / N}

  stopifnot(!any(diag(Cee) < 0))
  stopifnot(!any(diag(Cff) < 0))
  stopifnot(!any(diag(Cetet) < 0))
  stopifnot(!any(diag(Ceueu) < 0))
  return(
    list(
      X = X,
      Y = Y,
      Cxt = Cxt,
      Cxto = Cxto,
      Ctt = Ctt,
      Ctoto = Ctoto,
      Cyu = Cyu,
      Cyuo = Cyuo,
      Cuu = Cuu,
      Cuouo = Cuouo,
      Cee = Cee,
      Cff = Cff,
      Cetet = Cetet,
      Ceueu = Ceueu
    )
  )
}

M_step <- function(efit) {
  with(efit, {
    list(
      X = X,
      Y = Y,
      W = orth(Cxt %*% solve(Ctt))%*%diag(qr(Cxt)$qraux),
      C = orth(Cyu %*% solve(Cuu))%*%diag(qr(Cyu)$qraux),
      Wo = Cxto %*% solve(Ctoto),
      Co = Cyuo %*% solve(Cuouo),
      Phi1 = Cetet,
      Phi2 = Ceueu,
      Psi1 = Cee * diag(1,ncol(X)),
      Psi2 = Cff * diag(1,ncol(Y))
    )
  })
}

PJSC <- function(X, Y, n, nx, ny, nr_steps, tol){
  #set.seed(1)
  nx2 = max(nx,1)
  ny2 = max(ny,1)
  fito2m = o2m(X, Y, n, nx, ny, stripped=T)
  jhat=with(fito2m, list(u=W., v=C.))
  xhat=with(fito2m, list(v=P_Yosc.))
  yhat=with(fito2m, list(v=P_Xosc.))
  init = list(
    X = X,
    Y = Y,
    W = with(jhat,as.matrix(u)%*%diag(sign(colSums(as.matrix(u))),n)),
    C = with(jhat,as.matrix(v)%*%diag(sign(colSums(as.matrix(v))),n)),
    Wo = sign(nx)*with(xhat,as.matrix(v)%*%diag(sign(colSums(as.matrix(v))),nx2)),
    Co = sign(ny)*with(yhat,as.matrix(v)%*%diag(sign(colSums(as.matrix(v))),ny2)),
    Phi1 = diag(0.1, n),
    Phi2 = diag(0.1, n),
    Psi1 = diag(0.1 ^ 2, ncol(X)),
    Psi2 = diag(0.1 ^ 2, ncol(Y))
  )
  i = 0
  l_old = -Inf
  l_new = do.call(l_step, init)
  dvals = list(W = NULL, l = NULL)
  while (abs(l_new - l_old) > tol && i < nr_steps) {
    i = i + 1
    l_old = l_new
    fit = do.call(E_step, init)
    init = M_step(fit)
    print(crossprod(fit$Cxt %*% solve(fit$Ctt)))
    l_new = do.call(l_step, init)
    #dvals$W = c(dvals$W, cor(init$W[, 1], W[, 1]))
    dvals$l = c(dvals$l, l_old)
    if (l_new < l_old)
      message(paste("negative in step", i, ":", l_new - l_old))
  }
  message(paste("max=",which.max(dvals$l),"i=",i))
  list(
    Expec = fit,
    logl = dvals$l,
    est = init[-(1:2)]
  )
}

# heatmap.2(cor(cbind(X,Y)),Rowv=F,Colv=F, col=bluered(100))
# heatmap.2(cov2cor(do.call(cov_step, init)$joint),Rowv=F,Colv=F, col=bluered(100))
# heatmap.2(cov2cor(do.call(cov_step, init)$sys),Rowv=F,Colv=F, col=bluered(100))
# heatmap.2(cov2cor(cov(cbind(X,Y))-do.call(cov_step, init)$joint-do.call(cov_step, init)$sys),Rowv=F,Colv=F, col=bluered(100))
#
# #mapply(init[3:6], list(W=W,C=C,Wo=PYo,Co=PXo), FUN = cbind, SIMPLIFY = F)
# par(mfrow=c(3,2))
# invisible(sapply(mapply(list(W=W,C=C,Wo=PYo,Co=PXo),init[3:6], FUN = function(x,y) cbind(x[,1],y[,1]), SIMPLIFY = F), function(e) {plot(e);abline(a=0,b=1)}))
# with(dvals,{plot(W[-1]);plot(l[-1])})
# par(mfrow=c(1,1))

simulate_PJSC <- function(N, W, C, Wo, Co, Phi1, Phi2, Psi1, Psi2){
  r = ncol(W)
  rx = ncol(Wo)
  ry = ncol(Co)
  p = ncol(W)
  q = ncol(C)

  Z = matrix(rnorm(N*r))
  eps1 = mvrnorm(N*r, Sigma = Phi1, mu = rep(0,nrow(Phi1)))
  eps2 = mvrnorm(N*r, Sigma = Phi2, mu = rep(0,nrow(Phi2)))
  E = mvrnorm(N*p, Sigma = Psi1, mu = rep(0,nrow(Psi1)))
  Ff = mvrnorm(N*q, Sigma = Psi2, mu = rep(0,nrow(Psi2)))
  To = matrix(rnorm(N*rx))
  Uo = matrix(rnorm(N*ry))

  Tt = Z + eps1
  U = Z + eps2

  X = tcrossprod(Tt, W) + tcrossprod(To, Wo) + E
  Y = tcrossprod(U, C) + tcrossprod(Uo, Co) + Ff

  list(X = X, Y = Y,
       Tt = Tt, U = U, To = To, Uo = Uo,
       E = E, Ff = Ff, eps1 = eps1, eps2 = eps2)
}

