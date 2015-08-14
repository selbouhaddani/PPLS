#source('functions.R')
library(RcppEigen)
library(PPLS)
#dens= mvrnorm
#K=100
#EMsteps=500

set.seed(890)
N=1000
p=60;q=50
alp = 0.05

a = 2;
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

g = sigT[1]^2*B_T[1]^2 + sigH^2
Kw = sigT[1]^2 - sigT[1]^4*B_T[1]^2/sigY^2 + sigT[1]^4*B_T[1]^2*g/(sigY^2*(g+sigY^2))
Kc = g - sigT[1]^4*B_T[1]^2/sigX^2 + sigT[1]^6*B_T[1]^2/(sigX^2*(sigT[1]^2+sigX^2))
Kwc = sigT[1]^2*B_T[1]/(sigX^2*sigY^2) - Kc*sigT[1]^2*B_T[1]/(sigX^2*sigY^2*(Kc+sigY^2)) -
          sigT[1]^4*B_T[1]/(sigX^2*sigY^2*(sigT[1]^2+sigX^2)) +
          Kc*sigT[1]^4*B_T[1]/(sigX^2*sigY^2*(Kc+sigY^2)*(sigT[1]^2+sigX^2))
c1 = Kw / (sigX^2*(Kw + sigX^2)); c3 = Kc / (sigY^2*(Kc + sigY^2));

#Est = list()

#tic=proc.time();
#for(ind_i in 1:K){
Dat = simulC(N,W,C,PYo,PXo,B_T,sigX,sigY,sigH,sigT,0*sigTo,sigUo);
X = Dat[,1:p]
Y = Dat[,-(1:p)]

# Tt  = matrix(scale(rnorm(N*a)),N,a)%*%sigT
# TYo = matrix(scale(rnorm(N*a2)),N,a2)%*%sigTo
# UXo = matrix(scale(rnorm(N*b2)),N,b2)%*%sigUo
# #  E   = matrix(sigX*rnorm(N*p),N,p)
# E   = sigX*matrix(scale(rnorm(N*p)),N,p)
# #  Ff  = matrix(sigY*rnorm(N*q),N,q)#
# Ff  = sigY*matrix(scale(rnorm(N*q)),N,q)
# H_UT = sigH*matrix(scale(rnorm(N*a)),N,a)
#
# ## To model ###
# U = Tt %*% B_T
# U = U + H_UT
# X = tcrossprod(Tt,W) + E #+ tcrossprod(TYo,PYo)
# Y = tcrossprod(U,C) + Ff #+ tcrossprod(UXo,PXo)
# Dat = cbind(X,Y)

#ret = list(W=c(W),C=c(C),PYo=c(PYo),PXo=c(PXo),B=B_T,sig=c(sigX^2,sigY^2,sigH^2,sigT^2,sigTo^2,sigUo^2));
#ret = list(W=orth(runif(p)),C=orth(runif(q)),PYo=orth(runif(p)),PXo=orth(runif(q)),B=rchisq(1,1),sig=c(orth(rchisq(6,1))));
sim=o2m(X,Y,1,1,1)
ret = list(W=sim$W.,C=sim$C.,PYo=orth(sim$P_Yosc.),PXo=orth(sim$P_Xosc.),B=sim$B_T.,
           sig2=c(sum(diag(crossprod(sim$E)))/N/p,sum(diag(crossprod(sim$Ff)))/N/q,crossprod(sim$H_U)/N,
                 ssq(sim$Tt),ssq(sim$T_Y)*ssq(sim$P_Y),ssq(sim$U_X)*ssq(sim$P_X)/N))
ret2 = try(EMC(EMsteps,ret,Dat,0,0),T);
while(class(ret2)!="list"){
  #ret2=list(W=c(W*NA),C=c(C*NA),PYo=c(W*NA),PXo=c(C*NA))
  ret = list(W=orth(runif(p)),C=orth(runif(q)),PYo=orth(runif(p)),PXo=orth(runif(q)),B=rchisq(1,1),sig=c(orth(rchisq(6,1))));
  ret2 = try(EMC(EMsteps,ret,Dat,0,0),T);
}
Est[[ind_i]] = list(W=corr.o2m(ret2$W,W),C=corr.o2m(ret2$C,C),PYo=corr.o2m(ret2$PYo,PYo),PXo=corr.o2m(ret2$PXo,PXo))
#}
proc.time()-tic;
################################################################################

par(mfrow=c(2,2))
plot(ret2$W,ylab='W (X joint)',xlab='Variable');lines(W,col=2,lwd=2)
plot(ret2$PYo,ylab='P_Yo (X orth)',xlab='Variable');lines(PYo,col=2,lwd=2)
plot(ret2$C,ylab='C (Y joint)',xlab='Variable');lines(C,col=2,lwd=2)
plot(ret2$PXo,ylab='P_Xo (Y orth)',xlab='Variable');lines(PXo,col=2,lwd=2)
par(mfrow=c(1,1))

par(mfrow=c(2,2))
boxplot(t(sapply(1:K,function(i){Est[[i]]$W})));lines(W,col=2,lwd=2)
boxplot(t(sapply(1:K,function(i){Est[[i]]$PYo})));lines(PYo,col=2,lwd=2)
boxplot(t(sapply(1:K,function(i){Est[[i]]$C})));lines(C,col=2,lwd=2)
boxplot(t(sapply(1:K,function(i){Est[[i]]$PXo})));lines(PXo,col=2,lwd=2)
par(mfrow=c(1,1))
