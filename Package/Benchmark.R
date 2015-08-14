library(O2PLS)
library(microbenchmark)
tr <- function(X) sum(diag(X))
N = 600
p = 140
q = 35000
X = matrix(rnorm(N*p),N)
Y = matrix(rnorm(N*q),N)

library(proftools)
library(profr)
library(ggplot2)
library(lineprof)
Rprof_out <- Rprof(memory.profiling = TRUE);PPLS(X,Y,1,1e3,1e-4);Rprof(NULL)
plotProfileCallGraph(readProfileData(),score="total",nodeDetails=F,edgeDetails=T,edgesColored=T,layout="dot")


W = orth(rep(1,p))
C = orth(rep(1,q))

sigT = 1.05
B_T = 1.1
sigX = .05
sigY = .03
sigH = .01

g = sigT[1]^2*B_T[1]^2 + sigH^2
Kw = sigT[1]^2 - sigT[1]^4*B_T[1]^2/sigY^2 + sigT[1]^4*B_T[1]^2*g/(sigY^2*(g+sigY^2))
Kc = g - sigT[1]^4*B_T[1]^2/sigX^2 + sigT[1]^6*B_T[1]^2/(sigX^2*(sigT[1]^2+sigX^2))
Kwc = sigT[1]^2*B_T[1]/(sigX^2*sigY^2) - Kc*sigT[1]^2*B_T[1]/(sigX^2*sigY^2*(Kc+sigY^2)) -
  sigT[1]^4*B_T[1]/(sigX^2*sigY^2*(sigT[1]^2+sigX^2)) +
  Kc*sigT[1]^4*B_T[1]/(sigX^2*sigY^2*(Kc+sigY^2)*(sigT[1]^2+sigX^2))
c1 = Kw / (sigX^2*(Kw + sigX^2)); c3 = Kc / (sigY^2*(Kc + sigY^2));
c2 = Kwc

Sx = crossprod(X)
Sy = crossprod(Y)
Sxy = crossprod(X,Y)
Syx = crossprod(Y,X)
S = crossprod(cbind(X,Y))
all.equal(
tr(S%*%blockm(-c1*tcrossprod(W),-c2*tcrossprod(W,C),-c3*tcrossprod(C)) + S%*%diag(c(rep(1/sigX^2,p),rep(1/sigY^2,q))))
,
-(c1*ssq(X%*%W) + 2*c2*c(crossprod(X%*%W,Y%*%C)) + c3*ssq(Y%*%C)) + 1/sigX^2*ssq(X) + 1/sigY^2*ssq(Y)
)
