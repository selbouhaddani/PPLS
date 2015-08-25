func <- function(p,q,a,B_T,sigX,sigY,sigH,sigT){
  g = diag(sapply(1:a,function(i)sqrt(sigT[i,i]^2%*%B_T[i,i]^2 + sigH^2)),a)
  Kw = diag(sapply(1:a,function(i)sigT[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigY^2 + sigT[i,i]^4*B_T[i,i]^2*g[i,i]^2/(sigY^2*(g[i,i]^2+sigY^2))),a)
  Kc = diag(sapply(1:a,function(i)g[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigX^2 + sigT[i,i]^6*B_T[i,i]^2/(sigX^2*(sigT[i,i]^2+sigX^2))),a)
  Kwc = diag(sapply(1:a,function(i) sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2) - Kc[i,i]*sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)) -
                      sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(sigT[i,i]^2+sigX^2)) +
                      Kc[i,i]*sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)*(sigT[i,i]^2+sigX^2)) ),a)
  c1 = diag(sapply(1:a,function(i) Kw[i,i] / (sigX^2*(Kw[i,i] + sigX^2)) ),a)
  c3 = diag(sapply(1:a,function(i) Kc[i,i] / (sigY^2*(Kc[i,i] + sigY^2)) ),a)
  c2 = Kwc
  return(list(c1=c1,c2=c2,c3=c3,Kc=Kc))
}

f <- function(X,Y,W,C,B_T,sigX,sigY,sigH,sigT){
  p = nrow(W)
  q=nrow(C)
  a=ncol(W)
  
  g = diag(sapply(1:a,function(i)sqrt(sigT[i,i]^2%*%B_T[i,i]^2 + sigH^2)),a)
  Kw = diag(sapply(1:a,function(i)sigT[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigY^2 + sigT[i,i]^4*B_T[i,i]^2*g[i,i]^2/(sigY^2*(g[i,i]^2+sigY^2))),a)
  Kc = diag(sapply(1:a,function(i)g[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigX^2 + sigT[i,i]^6*B_T[i,i]^2/(sigX^2*(sigT[i,i]^2+sigX^2))),a)
  Kwc = diag(sapply(1:a,function(i) sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2) - Kc[i,i]*sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)) -
                      sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(sigT[i,i]^2+sigX^2)) +
                      Kc[i,i]*sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)*(sigT[i,i]^2+sigX^2)) ),a)
  c1 = diag(sapply(1:a,function(i) Kw[i,i] / (sigX^2*(Kw[i,i] + sigX^2)) ),a)
  c3 = diag(sapply(1:a,function(i) Kc[i,i] / (sigY^2*(Kc[i,i] + sigY^2)) ),a)
  c2 = Kwc
  
  logdet = sum(log(sigX^2+diag(sigT^2))) + (p-a)*log(sigX^2) + sum(log(sigY^2 + diag(Kc))) + (q-a)*log(sigY^2)
  TR = 1/sigX^2 * ssq(X) + 1/sigY^2 * ssq(Y)
  for(i in 1:a){TR = TR - c1[i,i] * ssq(X%*%W[,i]) - 2*c2[i,i]*crossprod(X%*%W[,i],Y%*%C[,i]) - c3[i,i]*ssq(Y%*%C[,i])}
  return(c(logdet,TR))
}

g = diag(sapply(1:a,function(i)sqrt(sigT[i,i]^2%*%B_T[i,i]^2 + sigH^2)),a)
Kw = diag(sapply(1:a,function(i)sigT[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigY^2 + sigT[i,i]^4*B_T[i,i]^2*g[i,i]^2/(sigY^2*(g[i,i]^2+sigY^2))),a)
Kc = diag(sapply(1:a,function(i)g[i,i]^2 - sigT[i,i]^4*B_T[i,i]^2/sigX^2 + sigT[i,i]^6*B_T[i,i]^2/(sigX^2*(sigT[i,i]^2+sigX^2))),a)
Kwc = diag(sapply(1:a,function(i) sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2) - Kc[i,i]*sigT[i,i]^2*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)) -
  sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(sigT[i,i]^2+sigX^2)) +
  Kc[i,i]*sigT[i,i]^4*B_T[i,i]/(sigX^2*sigY^2*(Kc[i,i]+sigY^2)*(sigT[i,i]^2+sigX^2)) ),a)
c1 = diag(sapply(1:a,function(i) Kw[i,i] / (sigX^2*(Kw[i,i] + sigX^2)) ),a)
c3 = diag(sapply(1:a,function(i) Kc[i,i] / (sigY^2*(Kc[i,i] + sigY^2)) ),a)
c2 = Kwc

A1 = diag(sigX^2,p) + tcrossprod(W%*%sigT)
A12 = tcrossprod(W%*%sigT^2,C%*%B_T)
A2 = diag(sigY^2,q) + tcrossprod(C%*%g)
A2inv = diag(1/sigY^2,q); for(i in 1:a) A2inv = A2inv - (g[i,i]^2)/(sigY^2*(g[i,i]^2+sigY^2))*tcrossprod(C[,i])
mse(A2inv, solve(A2))
A1inv = diag(1/sigX^2,p); for(i in 1:a) A1inv = A1inv - (sigT[i,i]^2)/(sigX^2*(sigT[i,i]^2+sigX^2))*tcrossprod(W[,i])
mse(A1inv, solve(A1))

C1 = diag(sigX^2,p); for(i in 1:a){C1 = C1 + Kw[i,i]*tcrossprod(W[,i])}
mse( C1 , A1 - A12 %*% A2inv %*% t(A12))

C2 = diag(sigY^2,q); for(i in 1:a){C2 = C2 + Kc[i,i]*tcrossprod(C[,i])}
mse( C2 , A2 - t(A12) %*% A1inv %*% A12)

mse(solve(sseXY_W()) , blockm(-W%*%c1%*%t(W)+diag(1/sigX^2,p) , -W%*%c2%*%t(C) , -C%*%c3%*%t(C)+diag(1/sigY^2,q)))




cl <- makeCluster(rep("localhost",detectCores()),type = "SOCK")
outp4 <- parSapply(cl,1:200,function(i){
  source('C:/Users/selbouhaddani/OneDrive/LUMC/PhD/Rcode/Projecten/PO2PLS/Rcpp/PreambleCpp.R')
  N=10000
  p=6;q=5
  alp = 0.05
  
  a = 1;
  a2 = b2 = a;
  
  #set.seed(2014)
  W   = orth(sapply(1:a,function(i){ dnorm(1:p,0.5*p+i/10*p,0.1*p) }),type=2)
  PYo = orth(sapply(1:a2,function(i){ dnorm(1:p,0.1*p+i/10*p,0.2*p) }))
  C   = orth(sapply(1:a,function(i){ dnorm(1:q,0.6*q+i/10*q,0.1*q) }),type=2)
  PXo = orth(sapply(1:b2,function(i){ dnorm(1:q,0.2*q+i/10*q,0.2*q) }))
  B_T = diag(exp(sapply(1:a,function(i)log(1.5)-(i-1)*.3)),a)#matrix(rnorm(a*b,0,25),nrow=a)
  
  sigT = diag(exp(sapply(1:a,function(i)0-(i-1)*.1)),a)
  sigTo = 0*diag(exp(sapply(1:a2,function(i)0-(i-1)*.1)),a2)
  sigUo = 0*diag(exp(sapply(1:b2,function(i)0-(i-1)*.1)),b2)
  sigX = .05#sqrt(alp/(1-alp)*(a*sigT^2+a2*sigTo^2)/p) #approx 5% var wrt X
  sigH = .05#sqrt(c(B_T^2*sigT^2))          #approx 5% var wrt U
  sigY = .05#%sqrt(alp/(1-alp)*(a*c(B_T^2)*sigT^2+a*sigH^2+b2*sigUo^2)/q) #approx 5% var wrt X
  
  #for(ind_i in 1:K){
  Dat = simulC(N,W,C,PYo,PXo,B_T,sigX,sigY,sigH,sigT,0*sigTo,0*sigUo);
  X = Dat[,1:p]
  Y = Dat[,-(1:p)]
  
  return(PPLS(X,Y,3,1e3)$Oth$Log)
  
})