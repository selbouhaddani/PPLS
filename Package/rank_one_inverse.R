f <- function(X,Y,W,C,B_T,sigX,sigY,sigH,sigT,c1,c2,c3,Kc){
  p = nrow(W)
  q=nrow(C)
  a=ncol(W)
  logdet = sum(log(sigX^2+diag(sigT^2))) + (p-a)*log(sigX^2) + sum(log(sigY^2 + diag(Kc))) + (q-a)*log(sigY^2)
  TR = 1/sigX^2 * ssq(X) + 1/sigY^2 * ssq(Y)
  for(i in 1:a){TR = TR - c1[i,i] * ssq(X%*%W[,i]) - 2*c2[i,i]*crossprod(X%*%W[,i],Y%*%C[,i]) - c3[i,i]*ssq(Y%*%C[,i])}
  return(TR)
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
