
##################################################################################
# This sets up necessary files to do EM on simulated data in Rcpp
##################################################################################
source('preambleCpp.R')

set.seed(890)
N=500
p=100;q=50
#p=500;q=250
alp = 0.5

a = 1;
a2 = b2 = 1;

dens= mvrnorm

K=100
EMsteps=500

#set.seed(2014)
W   = orth(sapply(1:a,function(i){ dnorm(1:p,0.5*p+i/10*p,0.1*p) }))
PYo = orth(sapply(1:a2,function(i){ dnorm(1:p,0.1*p+i/10*p,0.2*p) }))
C   = orth(sapply(1:a,function(i){ dnorm(1:q,0.6*q+i/10*q,0.1*q) }))
PXo = orth(sapply(1:b2,function(i){ dnorm(1:q,0.2*q+i/10*q,0.2*q) }))
B_T = diag(2,a)#matrix(rnorm(a*b,0,25),nrow=a)

sigT = 1
sigTo = 1
sigUo = 1
sigX = sqrt(alp/(1-alp)*(a*sigT^2+a2*sigTo^2)/p) #approx 5% var wrt X
sigH = sqrt(alp/(1-alp)*c(B_T^2*sigT^2))          #approx 5% var wrt U
sigY = sqrt(alp/(1-alp)*(a*c(B_T^2)*sigT^2+a*sigH+b2*sigUo^2)/q) #approx 5% var wrt X

Est = list()

tic=proc.time();
for(ind_i in 1:K){
Dat = simulC(N,W,C,PYo,PXo,c(B_T),sigX,sigY,sigH,sigT,sigTo,sigUo);
X = Dat[,1:p]
Y = Dat[,-(1:p)]

#ret = list(W=c(W),C=c(C),PYo=c(PYo),PXo=c(PXo),B=B_T,sig=c(sigX^2,sigY^2,sigH^2,sigT^2,sigTo^2,sigUo^2));
ret = list(W=orth(runif(p)),C=orth(runif(q)),PYo=orth(runif(p)),PXo=orth(runif(q)),B=rchisq(1,1),sig2=c(orth(rchisq(6,1))));

ret2 = try(EMC(EMsteps,ret,Dat,0.1,0.1),T);
while(class(ret2)!="list"){
  #ret2=list(W=c(W*NA),C=c(C*NA),PYo=c(W*NA),PXo=c(C*NA))
  ret = list(W=orth(runif(p)),C=orth(runif(q)),PYo=orth(runif(p)),PXo=orth(runif(q)),B=rchisq(1,1),sig2=c(orth(rchisq(6,1))));
  ret2 = try(EMC(EMsteps,ret,Dat,0.1,0.1),T);
}
Est[[ind_i]] = list(W=corr.o2m(ret2$W,W),C=corr.o2m(ret2$C,C),PYo=corr.o2m(ret2$PYo,PYo),PXo=corr.o2m(ret2$PXo,PXo))
}
proc.time()-tic;
################################################################################

par(mfrow=c(2,2))
plot(ret2$W,ylab='W (X joint)',xlab='Variable');lines(W,col=2,lwd=2)
plot(-ret2$PYo,ylab='P_Yo (X orth)',xlab='Variable');lines(PYo,col=2,lwd=2)
plot(ret2$C,ylab='C (Y joint)',xlab='Variable');lines(C,col=2,lwd=2)
plot(ret2$PXo,ylab='P_Xo (Y orth)',xlab='Variable');lines(PXo,col=2,lwd=2)
par(mfrow=c(1,1))