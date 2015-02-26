multiroot2 = function(f, start, Cxt,Cxto,Ctt,Ctto,Ctoto)
{
  outp = multiroot(f, start,Cxt=Cxt,Cxto=Cxto,Ctt=Ctt,Ctto=Ctto,Ctoto=Ctoto);
  sol = 12345;
  if(is.finite(outp$estim) && abs(outp$estim) < 1e-4){sol = outp$root};
  return(sol);
}

sseXY <- function()
  {
  SX. = tcrossprod(W)+tcrossprod(PYo)+sigX^2*diag(1,nrow=p)
  SY. = tcrossprod(C%*%t(B_T))+tcrossprod(C)*sigH^2+tcrossprod(PXo)+sigY^2*diag(1,nrow=q)
  SXY. = W%*%B_T%*%t(C)
  #invSX. = woodinv(sigX.,cbind(W.,P_Yosc.))
  #invSY. = woodinv(sigY.,cbind(C.,P_Xosc.))
  Sfull = blockm(SX.,SXY.,SY.)
  return(Sfull)
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

corr.o2m = function(W_hat,W){
  sign_What = sign(c(crossprod(W_hat,W)))
  return(sign_What*W_hat)
}

woodinv <- function(sig,Uw)
  #input: sigma of noise, loading U
  #output: inverse of (U %*% t(U) + sig*Ip) via Woodbury identity
{
  p = length(Uw[,1])
  p2 = length(Uw[1,])
  Ai = diag(1/sig,p)
  invS = Ai - Ai%*%Uw%*%solve(diag(1,p2)+t(Uw)%*%Ai%*%Uw)%*%t(Uw)%*%Ai
  return(invS)
}

logl <- function(W.=W, C.=C, P_Yosc.=PYo, P_Xosc.=PXo, B_T.=B_T, X.=X, Y.=Y,
                 sigX.=sigX,sigY.=sigY,sigH.=sigH,sigT.=sigT,sigTo.=sigTo,sigUo.=sigUo)
  #input: Loadings, parameters
  #output: log-likelihood of Probabilistic O2PLS
{
  Data=cbind(X.,Y.)
  N = dim(X.)[1]
  p = dim(X.)[2]
  q = dim(Y.)[2]
  
  # Sigma X
  SX. = sigT.^2* tcrossprod(W.)+sigTo.^2*tcrossprod(P_Yosc.)+sigX.^2*diag(1,nrow=p)
  SXY. = sigT.^2* W.%*%B_T.%*%t(C.)
  SY. = sigT.^2* tcrossprod(C.%*%t(B_T.))+tcrossprod(C.)*sigH.^2+sigUo.^2*tcrossprod(P_Xosc.)+sigY.^2*diag(1,nrow=q)
  Sobs. = rbind( cbind(SX.,SXY.) , cbind(t(SXY.),SY.) )
  
  invSobs. = solve(Sobs.)
  #som = 0
  #for(i in 1:N){som = som + t(Data[i,])%*%invSobs.%*%t(t(Data[i,]))}
  #som=sum(sapply(1:N,function(i){t(Data[i,])%*%invSobs.%*%t(t(Data[i,]))}))
  som = sum(diag(crossprod(cbind(X.,Y.))%*%invSobs.))
  logdet = c(determinant(Sobs.)$mod)
  Q = - N*0.5*logdet - 0.5*som
  return(Q)
}

maxl = function(W)
  #input: Vector W
  #output: a +/- 1 on the absolute maximal element of W, zeros elsewhere
{
  return(sign(W)*(abs(W)==max(abs(W))))
}

regr<-function(X,Y)
  #input: prediction matrix T (Nxa) and response matrix X(Nxp)
  #output: Yields loadings/coefs for X~T (pxa)
{
  return(crossprod(Y,X)%*%solve(crossprod(X)))
}

orth<-function(X)
{
  #input is some matrix X
  #output is a matrix with orthonormal columns (via SVD(X))
  e=svd(X)
  return(tcrossprod(e$u,e$v))
}

# b.solve<-function(A,B,C,D){
#   # invert block matrices [A,B \\ C,D]
#   Inv = solve(rbind(cbind(A,B),cbind(C,D)))
#   return(Inv)
# }

mse <- function(x,y)
  #Input: 2 matrices x,y
  #Output: mean squared distance between 2 matrices (number)
{
  mean((x-y)^2)
}

boxtr<-function(X,L=0)
  #Input: matrix X, box-cox parameter L
  #Output: Box-Cox transform of X
{
  if(L==0){return(log(X))}
  if(L!=0){return((X^L-1)/L)}
}

pca <- function(X,ncomp=3,scale_var=T)
  #Input: Matrix X; number of PC's wanted; scale?
  #Output: PCA(X) with scores; loadings; relative variation of each PC
{
  Xpca=X
  if(scale_var){Xpca=scale(X,scale=F)}   #Scale
  Xsvd=svd(Xpca,nu=ncomp,nv=ncomp)       #perform SVD
  
  #distinguish whether X is a vector or matrix
  #Defines T = UD; P = V
  Xsvd.sc=NA
  if(ncomp>0){
    if(is.matrix(X)&&!any(dim(as.matrix(X))==1)){Xsvd.sc=Xsvd$u %*% diag(Xsvd$d[1:ncomp],ncomp)}
    if(!is.matrix(X)||any(dim(as.matrix(X))==1)){Xsvd.sc=Xsvd$u * Xsvd$d}
  }
  # D^2 contains variances of each PC
  vars=Xsvd$d^2 
  
  #List scores loadings and relative variances
  list(scores=Xsvd.sc,loadings=Xsvd$v,ev=Xsvd$d,variation=cumsum(vars/sum(vars)))
}

pcaplot<-function(X,richt=1:2,amp=1,plotmax=F,load=NULL,...)
  #Input: Datamatrix X; which X-variable (1st, 2nd,); amplifier of arrow for visibility
  #       plot direction of max loading?; loading matrix
  #Output: Plot of (X_i,X_j) with arrow pointing in wanted direction
{
  # potential X_i and X_j
	i=richt[1];j=richt[2]
  
  # Calculate loadings from svd or using given load
  if(is.null(load)){P=svd(X,nu=0,nv=1)$v}
	if(!is.null(load)){P=load}
  
  # Calculates X-variable with max loading if plotmax=T
  if(plotmax){i=which.max(abs(P)[,1]);j=which.max(abs(P)[-i,1])}
  
  #Define X_i and X_j and arrows and plot
	X1=X[,i]
	X2=X[,j]
  vec1=cumsum(c(mean(X1),amp*P[i,1]))
  vec2=cumsum(c(mean(X2),amp*P[j,1]))
	plot(X[,c(i,j)],xlab=paste("var ",i),ylab=paste("var ",j),...)
	arrows(x0=vec1[1],y0=vec2[1],x1=vec1[2],y1=vec2[2],col=2,lwd=2)
  
  #vec1 is x-coord of endpoint of arrow, vec2 is y-coord
  list(vec1=vec1,vec2=vec2)
}

plsm <- function(X,Y,ncomp=3,scale_var=F)
  #Input: Data X,Y; number of wanted components; scale in variance?
  #Output: class 'plsm' list with statistics, regression with intercept
  {
  meanX=colMeans(X)
  W=NULL;Tt=NULL;P=NULL;Q=NULL;
  E=X#scale(X,scale=scale_var);          #scale X
  Ff=Y#scale(Y,scale=scale_var)          #scale Y
  totvar=ssq(E)                          #total variation of scaled X
  if(ncomp>min(dim(X))){stop(paste("Check: ncomp <=",min(dim(X))))}
  
  #Extract PLS components
  for(i in 1:ncomp){
    svdec=svd(crossprod(E,Ff),nu=1,nv=0)#XtY=UDV
    w=svdec$u[,1]                      #1st PC for XtY
    t=E%*%w                            #corresponding score vector
    p=crossprod(E,t)/c(crossprod(t))   #regressing E on t
    q=crossprod(Ff,t)/c(crossprod(t))  #regressing F on t
    E=E-tcrossprod(t,p)                #remove current PC
    Ff=Ff-tcrossprod(t,q)              #remove current direction
    W=cbind(W,w);Tt=cbind(Tt,t)        #bind all vectors
    P=cbind(P,p);Q=cbind(Q,q)          #bind all vectors
  }
  
  #Further relevant matrices/coefficients
  A=solve(crossprod(cbind(1,Tt)))%*%crossprod(cbind(1,Tt),Y) #coef incl intercept
  R=W%*%solve(crossprod(P,W))          #Generalized P^-1
  B=R%*%(as.matrix(A)[-1,])            #coef in original space excl intercept
  interc=apply(as.matrix(Y),2,mean)-apply(as.matrix(X%*%B),2,mean)
  B=rbind(interc,B)
  
  #per PC total variation; Ytotal variation; variation in Y explained by X
  vars=c(sum(diag(crossprod(P[,1])*crossprod(Tt[,1]))),diff(cumsum(diag(crossprod(P)*crossprod(Tt)))))
  Ytotvar=diag(crossprod(scale(Y,scale=scale_var)))
  Yvars=diag(crossprod(Tt[,1]%*%t(as.matrix(A)[2,])))/Ytotvar
  if(ncomp>1){
    for(j in 2:ncomp){
      Yvars=rbind( Yvars,diag(crossprod( Tt[,1:j]%*%as.matrix(A)[2:(j+1),] ))/Ytotvar )
    }
  }
  
  fit.val=cbind(1,X)%*%B
  model=list(scores=Tt,Xmean=meanX,loadings=P,Pinv=R,coef=A,origcoef=B,Xweights=W,Yweights=Q,
	  variation=cumsum(vars)/totvar,Yvariation=Yvars,fitted.values=fit.val)
  class(model)='plsm'
  return(model)
}

ssq<-function(X)
  #Input: matrix X
  #Output: sum of squares of X
{
  return(sum(X^2))
}

o2m<-function(X,Y,n,nx,ny)
  #Input: Data X,Y; number of correlated and orthogonal components (nx is dim(X_orth))
  #Output: O2PLS class 'o2m' model with statistics, (Trygg, O2PLS)
{
  X_true = X
  Y_true = Y
  
  T_Yosc = U_Xosc = matrix(0,length(X[,1]),1)
  P_Yosc = W_Yosc = matrix(0,length(X[1,]),1)
  P_Xosc = C_Xosc = matrix(0,length(Y[1,]),1)
  
  Y = as.matrix(Y)
  N=dim(X)[1]
  p=dim(X)[2]; q=dim(Y)[2]
  
  if(nx*ny>0){
    # number of latent components is n-max(nx,ny), so make sure n>max(nx,ny)
    n2=n+max(nx,ny)
    
    cdw = svd(t(Y)%*%X,nu=n2,nv=n2); # 3.2.1. 1
    C=cdw$u;W=cdw$v
    
    Tt = X%*%W;                    # 3.2.1. 2
    
    if(nx > 0){
      # 3.2.1. 3  Note here the E_XY should be X - TW' 
      E_XY = X - Tt%*%t(W);
      
      # remove orthogonal components from X sequentially
      # for lv = 1:LVX
      #     [w_Yosc,s,v] = svds(E_XY'*T,1);
      #     t_Yosc = X*w_Yosc;
      #     p_Yosc = (inv(t_Yosc'*t_Yosc)*t_Yosc'*X)';
      #     X = X - t_Yosc*p_Yosc';
      #     
      #     # Collect the Y orthogonal X components
      #     T_Yosc = [T_Yosc t_Yosc];
      #     P_Yosc = [P_Yosc p_Yosc];
      #     W_Yosc = [W_Yosc w_Yosc];
      # end
      # Even niet sequentieel
    
      udv = svd(t(E_XY)%*%Tt,nu=nx,nv=0);
      W_Yosc = udv$u ; s = udv$d 
      T_Yosc = X%*%W_Yosc;
      P_Yosc = t(solve(t(T_Yosc)%*%T_Yosc)%*%t(T_Yosc)%*%X);
      X = X - T_Yosc%*%t(P_Yosc);
      
      # Update T again (since X has changed)
      Tt = X%*%W;
    }
    
    U = Y%*%C;                   # 3.2.1. 4
    
    if(ny > 0){
      # 3.2.1. 5
      F_XY = Y - U%*%t(C);
                      
      # calculate LVY orthogonal Y components
      # for lv = 1:LVY
      #     [c_Xosc,s,v] = svds(F_XY'*U,1);
      #     u_Xosc = Y*c_Xosc;
      #     p_Xosc = (inv(u_Xosc'*u_Xosc)*u_Xosc'*Y)';
                      #     Y = Y - u_Xosc*p_Xosc';   
      #     
      #     U_Xosc = [U_Xosc u_Xosc];
      #     P_Xosc = [P_Xosc p_Xosc];
      #     C_Xosc = [C_Xosc c_Xosc];
      # end 
      # Even niet sequentieel
      
      udv = svd(t(F_XY)%*%U,nu=ny,nv=0);
      C_Xosc = udv$u ; s = udv$d 
      U_Xosc = Y%*%C_Xosc;
      P_Xosc = t(solve(t(U_Xosc)%*%U_Xosc)%*%t(U_Xosc)%*%Y);
      Y = Y - U_Xosc%*%t(P_Xosc);
    
      # Update U again (since Y has changed) 
      U = Y%*%C;
    }
  }
  # repeat steps 1, 2, and 4 before step 6
  cdw = svd(t(Y)%*%X,nu=n,nv=n);    # 3.2.1. 1
  C=cdw$u;W=cdw$v
  Tt = X%*%W;                    # 3.2.1. 2
  U = Y%*%C;                    # 3.2.1. 4
  
  # 3.2.1. 6
  B_U = solve(t(U)%*%U)%*%t(U)%*%Tt;
  B_T = solve(t(Tt)%*%Tt)%*%t(Tt)%*%U;
  
  # Other model components
  E = X_true - Tt%*%t(W) - T_Yosc%*%t(P_Yosc);
  Ff = Y_true - U%*%t(C) - U_Xosc%*%t(P_Xosc);
  H_TU = Tt - U%*%B_U;
  H_UT = U - Tt%*%B_T;
  Y_hat = Tt%*%B_T%*%t(C);
  X_hat = U%*%B_U%*%t(W);
  
  R2Xcorr = (ssq(Tt%*%t(W))/ssq(X_true))
  R2Ycorr = (ssq(U%*%t(C))/ssq(Y_true))
  R2X_YO = (ssq(T_Yosc%*%t(P_Yosc))/ssq(X_true))
  R2Y_XO = (ssq(U_Xosc%*%t(P_Xosc))/ssq(Y_true))
  R2Xhat = 1 - (ssq(U%*%B_U%*%t(W) - X_true)/ssq(X_true))
  R2Yhat = 1 - (ssq(Tt%*%B_T%*%t(C) - Y_true)/ssq(Y_true))
  R2X = R2Xcorr + R2X_YO
  R2Y = R2Ycorr + R2Y_XO
  
  model=list(
    Tt=Tt,W.=W,U=U,C.=C,E=E,Ff=Ff,T_Yosc=T_Yosc,P_Yosc.=P_Yosc,W_Yosc=W_Yosc,U_Xosc=U_Xosc,P_Xosc.=P_Xosc,
    C_Xosc=C_Xosc,B_U=B_U,B_T.=B_T,H_TU=H_TU,H_UT=H_UT,X_hat=X_hat,Y_hat=Y_hat,R2X=R2X,R2Y=R2Y,
    R2Xcorr=R2Xcorr,R2Ycorr=R2Ycorr,R2X_YO=R2X_YO,R2Y_XO=R2Y_XO,R2Xhat=R2Xhat,R2Yhat=R2Yhat
  )
  class(model)='o2m'
  return(model)
}

o2mV2<-function(X,Y,n,nx,ny)
  #Input: Data X,Y; number of correlated and orthogonal components (nx is dim(X_orth))
  #Output: O2PLS class 'o2m' model with statistics, (Trygg, O2PLS)
  # Unique model components!
{
  Y = as.matrix(Y)
  N=dim(X)[1]
  p=dim(X)[2]; q=dim(Y)[2]
  
  T_Yosc = U_Xosc = matrix(0,length(X[,1]),1)
  P_Yosc = W_Yosc = matrix(0,length(X[1,]),1)
  P_Xosc = C_Xosc = matrix(0,length(Y[1,]),1)
  
  X_true = X
  Y_true = Y
  
  cdw = svd(t(Y)%*%X,nu=n+max(nx,ny),nv=n+max(nx,ny)); # 3.2.1. 1
  if(length(Y[1,])<length(X[1,])){
    C2 = cdw$u;W=cdw$v
    C = C2 %*% diag(sign(C2[1,]),n+max(nx,ny)) 
    W = W %*% diag(sign(C2[1,]),n+max(nx,ny))
  }# make decomp unique in minus sign
  if(length(Y[1,])>=length(X[1,])){
    C = cdw$u;W2=cdw$v
    C = C %*% diag(sign(W2[1,]),n+max(nx,ny)) 
    W = W2 %*% diag(sign(W2[1,]),n+max(nx,ny))
  }
  Tt = X%*%W;                    # 3.2.1. 2
  
  # 3.2.1. 3  Note here the E_XY should be X - TW' 
  E_XY = X - Tt%*%t(W);
  
  # remove orthogonal components from X sequentially
  # for lv = 1:LVX
  #     [w_Yosc,s,v] = svds(E_XY'*T,1);
  #     t_Yosc = X*w_Yosc;
  #     p_Yosc = (inv(t_Yosc'*t_Yosc)*t_Yosc'*X)';
  #     X = X - t_Yosc*p_Yosc';
  #     
  #     # Collect the Y orthogonal X components
  #     T_Yosc = [T_Yosc t_Yosc];
  #     P_Yosc = [P_Yosc p_Yosc];
  #     W_Yosc = [W_Yosc w_Yosc];
  # end
  # Even niet sequentieel
  if(nx > 0){
    udv = svd(t(E_XY)%*%Tt,nu=nx,nv=0);
    W_Yosc = udv$u ; s = udv$d 
    W_Yosc = W_Yosc %*% diag(sign(W_Yosc[1,]),nx)
    T_Yosc = X%*%W_Yosc;
    P_Yosc = t(solve(t(T_Yosc)%*%T_Yosc)%*%t(T_Yosc)%*%X);
    X = X - T_Yosc%*%t(P_Yosc);
    
    svdP = svd(P_Yosc)
    Cosc = svdP$v%*%diag(svdP$d,length(svdP$d))%*%t(svdP$v)
    Cosc_inv = svdP$v%*%diag(1/svdP$d,length(svdP$d))%*%t(svdP$v)
    P_Yosc = t(Cosc_inv%*%t(P_Yosc))
    T_Yosc = T_Yosc%*%Cosc
  }
  # Update T again (since X has changed)
  Tt = X%*%W;
  
  U = Y%*%C;                   # 3.2.1. 4
  
  # 3.2.1. 5
  F_XY = Y - U%*%t(C);
  
  # calculate LVY orthogonal Y components
  # for lv = 1:LVY
  #     [c_Xosc,s,v] = svds(F_XY'*U,1);
  #     u_Xosc = Y*c_Xosc;
  #     p_Xosc = (inv(u_Xosc'*u_Xosc)*u_Xosc'*Y)';
  #     Y = Y - u_Xosc*p_Xosc';   
  #     
  #     U_Xosc = [U_Xosc u_Xosc];
  #     P_Xosc = [P_Xosc p_Xosc];
  #     C_Xosc = [C_Xosc c_Xosc];
  # end 
  # Even niet sequentieel
  if(ny > 0){
    udv2 = svd(t(F_XY)%*%U,nu=ny,nv=0);
    C_Xosc = udv2$u ; s = udv2$d 
    C_Xosc = C_Xosc %*% diag(sign(C_Xosc[1,]),ny)
    U_Xosc = Y%*%C_Xosc;
    P_Xosc = t(solve(t(U_Xosc)%*%U_Xosc)%*%t(U_Xosc)%*%Y);
    Y = Y - U_Xosc%*%t(P_Xosc);
    
    svdP = svd(P_Xosc)
    Cosc = svdP$v%*%diag(svdP$d,length(svdP$d))%*%t(svdP$v)
    Cosc_inv = svdP$v%*%diag(1/svdP$d,length(svdP$d))%*%t(svdP$v)
    P_Xosc = t(Cosc_inv%*%t(P_Xosc))
    U_Xosc = U_Xosc%*%Cosc
  }
  
  # Update U again (since Y has changed) 
  U = Y%*%C;
  
  # repeat steps 1, 2, and 4 before step 6
  udv3 = svd(t(Y)%*%X,nu=n,nv=n);    # 3.2.1. 1
  if(length(Y[1,])<length(X[1,])){
  C2 = udv3$u;W=udv3$v
  C = C2 %*% diag(sign(C2[1,]),n) 
  W = W %*% diag(sign(C2[1,]),n)
  }# make decomp unique in minus sign
  if(length(Y[1,])>=length(X[1,])){
  C = udv3$u;W2=udv3$v
  C = C %*% diag(sign(W2[1,]),n) 
  W = W2 %*% diag(sign(W2[1,]),n)
  }
  Tt = X%*%W;                    # 3.2.1. 2
  U = Y%*%C;                    # 3.2.1. 4
  
  # 3.2.1. 6
  B_U = solve(t(U)%*%U)%*%t(U)%*%Tt;
  B_T = solve(t(Tt)%*%Tt)%*%t(Tt)%*%U;
  
  # Other model components
  E = X_true - Tt%*%t(W) - T_Yosc%*%t(P_Yosc);
  Ff = Y_true - U%*%t(C) - U_Xosc%*%t(P_Xosc);
  H_TU = Tt - U%*%B_U;
  H_UT = U - Tt%*%B_T;
  Y_hat = Tt%*%B_T%*%t(C);
  X_hat = U%*%B_U%*%t(W);
  
  R2Xcorr = (ssq(Tt%*%t(W))/ssq(X_true))
  R2Ycorr = (ssq(U%*%t(C))/ssq(Y_true))
  R2X_YO = (ssq(T_Yosc%*%t(P_Yosc))/ssq(X_true))
  R2Y_XO = (ssq(U_Xosc%*%t(P_Xosc))/ssq(Y_true))
  R2Xhat = 1 - (ssq(U%*%B_U%*%t(W) - X_true)/ssq(X_true))
  R2Yhat = 1 - (ssq(Tt%*%B_T%*%t(C) - Y_true)/ssq(Y_true))
  R2X = R2Xcorr + R2X_YO
  R2Y = R2Ycorr + R2Y_XO
  
  model=list(
    Tt=Tt,W.=W,U=U,C.=C,E=E,Ff=Ff,T_Yosc=T_Yosc,P_Yosc.=P_Yosc,W_Yosc=W_Yosc,U_Xosc=U_Xosc,P_Xosc.=P_Xosc,
    C_Xosc=C_Xosc,B_U=B_U,B_T.=B_T,H_TU=H_TU,H_UT=H_UT,X_hat=X_hat,Y_hat=Y_hat,R2X=R2X,R2Y=R2Y,
    R2Xcorr=R2Xcorr,R2Ycorr=R2Ycorr,R2X_YO=R2X_YO,R2Y_XO=R2Y_XO,R2Xhat=R2Xhat,R2Yhat=R2Yhat
  )
  class(model)='o2m'
  return(model)
}

summary.o2m<-function(fit)
  #Input: O2PLS model (class o2m)
  #Output: table/matrix with R2 statistics
  {
  a=length(fit$W.[1,])
  Mname=list(c(""),c("Comp","R2X","R2Y","R2Xcorr",'R2Ycorr','R2Xhat',"R2Yhat","XRatio","YRatio"))
  M=matrix(c(a/100,fit$R2X,fit$R2Y,fit$R2Xcorr,fit$R2Ycorr,fit$R2Xhat,fit$R2Yhat,fit$R2Xhat/fit$R2Xcorr,fit$R2Yhat/fit$R2Ycorr),nrow=1,dimnames=Mname)
  return(round(100*M,2))
}

rmsep <- function(Xtst,Ytst,fit,n_orth=1)
  #Input: test matrices (now: row i of data X,Y); used model
  #Output: MSE of Yhat from model, wrt Ytest
{
  n1=n2=T
  if(!is.matrix(Xtst)){Xtst=t(Xtst);n1=F} # If Xtst is a row-vector
  if(!is.matrix(Ytst)){Ytst=t(Ytst);n2=F} # If Xtst is a row-vector

  if(class(fit)=='plsm')
    # PLS fit is quickly done
  {
    Yhat=cbind(1,Xtst)%*%fit$ori
  }
  
  if(class(fit)=='oplsm')
  # OPLS corrects Xtst
  {
    if(n1){Xtst=oscr(Xtst,Ytst,n_orth=n_orth)$Xcorr}
    Yhat=cbind(1,Xtst)%*%fit$ori
  }
  
  if(class(fit)=='o2m')
  # O2PLS we should correct both Xtst and Ytst?
  {
      if(n1)
    {
      #To=Xtst%*%fit$W_Yosc
      #Po=t(solve(t(To)%*%To)%*%t(To)%*%Xtst)
      Xtst=Xtst# - tcrossprod(To,Po)
    }
      Yhat = Xtst%*%fit$W.%*%fit$B_T%*%t(fit$C.)
      Xhat = Ytst%*%fit$C.%*%fit$B_U%*%t(fit$W.)
  }
  return(mean(c(sqrt(mse(Yhat,Ytst)))))#,sqrt(mse(Xhat,Xtst)))))
}

mserr = function(i,a,Datalist,a2,b2)
  #Input: i-th row to exclude, a, list or data.frame with X and Y, nx, ny
  #Output: rmse of Yhat_i and Y[i]
{
  X = Datalist$X
  Y = Datalist$Y
  pars = list(X = X[-i,] , Y = Y[-i,] , n = a , nx = a2 , ny = b2)
  fit=try(do.call(o2m,pars) , silent=T)
  return( ifelse(class(fit)=="try-error" , NA , rmsep(X[i,],Y[i,],fit)) )
}

loocv <- function(X,Y,a=1:2,a2=1,b2=1,fitted_model=NULL,func=o2m_stripped,app_err=F,kcv)
  #Input: Data X,Y; model to fit; loop through nr of components; 
          # calculate apparent error?; nr of folds (loo:kcv=N)
  #Output: several MSE's per chosen nr of component
{
  if(!is.null(fitted_model)){app_err=F;warning("apparent error calculated with provided fit")}
  p=dim(X)[2]
  # determine type of model
  type=3#ifelse(deparse(substitute(func))=="o2m",3,ifelse(deparse(substitute(func))=="oplsm",2,1))
  
  N = length(X[,1]);if(N!=length(Y[,1])){stop('N not the same')}
  mean_err=mean_fit=NA*1:max(length(a),length(a2),length(b2))
  k=0
  
  #blocks contains the begin and endpoints of test indices to use
  blocks = c(seq(0,N,by=floor(N/kcv)),N)
  
  #loop through chosen parameters
  for(j in a){for(j2 in a2){for(j3 in b2){
  k=k+1
  err=NA*1:kcv
  folds=sample(N)
  #loop through number of folds
  for(i in 1:kcv){
    ii = (blocks[i]+1):(blocks[i+1])
    if(type==3){pars=list(X=X[-folds[ii],],Y=Y[-folds[ii],],n=j,nx=j2,ny=j3)}
    if(type==2){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j,n_orth=j2)}
    if(type==1){pars=list(X=X[-i,],Y=Y[-i,],ncomp=j)}
      fit=try(do.call(func,pars),silent=T)
      err[i] = ifelse(class(fit)=="try-error",NA,rmsep(X[folds[ii],],Y[folds[ii],],fit))
  }
  mean_err[k]=mean(err)
  #calculate apparent error
  if(app_err && is.null(fitted_model)){
    if(class(fit)=='o2m'){pars2=list(X=X,Y=Y,n=j,nx=j2,ny=j3)}
    if(class(fit)=='oplsm'){pars2=list(X=X,Y=Y,ncomp=j,n_orth=j2)}
    if(class(fit)=='plsm'){pars2=list(X=X,Y=Y,ncomp=j)}
    fit2=try(do.call(func,pars2),T)
    mean_fit[k]=ifelse(class(fit)=="try-error",NA,rmsep(X,Y,fit2))
    print('1e loop')
  }
  if(!is.null(fitted_model))
  {
    mean_fit[k]=rmsep(X,Y,fitted_model)
    print('2e loop')
  }
}}}
return(list(CVerr=mean_err,Fiterr=mean_fit))
}

adjR2 <- function(X,Y,a=1:2,a2=1,b2=1,parall=F,cl=NULL)
  #Input: Data X,Y; model to fit; loop through nr of components; 
  # calculate apparent error?; nr of folds (loo:kcv=N)
  #Output: several MSE's per chosen nr of component
{
  X;Y;
  cl_was_null=FALSE
  if(!parall){S_apply=function(cl=NULL,x,fun){sapply(x,fun)}}
  if(parall&is.null(cl)){
    cl_was_null=TRUE
    S_apply=parSapply
    cl <- makeCluster(rep('localhost', detectCores()),type='SOCK')
    clusterExport(cl=cl, varlist=c("ssq",'o2m_stripped',"adjR2"))
  }
  if(parall&!is.null(cl)){S_apply=parSapply}
  
  pars1 = merge(merge(data.frame(a = a),data.frame(a2=a2)),data.frame(b2=b2))
  pars2 = apply(pars1, 1, as.list)
  N = dim(X)[1]
  #cl <- makeCluster(rep( 'localhost', detectCores()),type='SOCK')
  #clusterExport(cl=cl, varlist=c("X","Y","N","pars2","ssq",'o2m'))
  outp=S_apply(cl,pars2,function(p){
    fit=o2m_stripped(X,Y,p$a,p$a2,p$b2);
    RX = 1-ssq(fit$H_UT)/ssq(fit$U)
    RY = 1-ssq(fit$H_TU)/ssq(fit$Tt)
    adjRX = RX#1 - (1 - RX)*(N - 1)/(N - p$a - 1)
    adjRY = RY#1 - (1 - RY)*(N - 1)/(N - p$a - 1)
    return(list(adjR2X = adjRX, adjR2Y = adjRY))
    })
  if(parall&cl_was_null==TRUE){stopCluster(cl)}
  return(outp)
}

vnorm <- function(x)
  #Input: matrix x
  #Output: 2-norm of X per column ( ssq(X) = sum(vnorm(X)) )
{
  x = as.matrix(x)
  return(sqrt(apply(x^2,2,sum)))
}

BCquant=function(boot,est)
  #Input: bootstrap samples matrix p x B, estimate vector p x 1  
  #Output: list of lower alpha's and upper alpha's 
  #        (usually different from 0.025,0.975)
{
  est=as.matrix(est)
  boot=as.matrix(boot)
  p = length(boot[,1])
  z0 = sapply(1:p,function(i){qnorm(mean(est[i,]>boot[i,]))})
  low=sapply(1:p,function(i){pnorm(2*z0[i]-1.96)})
  up=sapply(1:p,function(i){pnorm(2*z0[i]+1.96)})
  quant=rbind(lower=low,upper=up)
  return(sapply(1:p,function(i){quantile(boot[i,],quant[,i])}))
}

momentsoutlier<-function(X)
  #Input: Data X
  #Output: skewness and kurtosis analysis on X
{
  # n is number of rows
  n=length(X[,1])
  
  #SE. is Standard error of measure (Skewness or Kurtosis)
  SES=(6*n*(n-1))/((n-2)*(n+1)*(n+3))
  SES=sqrt(SES)
  
  #.ind is indices of X-variables that are > 3*SE
  Sind=which(abs(apply(X,2,skewness))>3*SES)
  
  SEK=(n^2-1)/((n-3)*(n+5))
  SEK=sqrt(SEK)
  Kind=which(abs(apply(X,2,kurtosis))>3*SEK)
  
  list(SES=SES,SEK=SEK,n=n,S.ind=Sind,K.ind=Kind)
}

momentsplot<-function(X,type=1:4,ident=F,corr=F,SES=F,SEK=SES)
  #Input: Data X, type of plot (location, spread, skewness, kurtosis)
          # identify manually outliers? plot normal X vs log(X) if type==3?
  #Output: plots with the 4 characteristics
{
  #initiate plot par
  l=length(type)
  if(l<4){par(mfrow=c(l,1))}
  if(l==4){par(mfrow=c(l-2,2))}
  
  #identify only makes sense if 1 plot is asked
  #ident=(l==1)*ident
  
  #Plot mean - median
  if(sum(type==1)>0){
    X2=scale(X)
    plot(apply(X2,2,mean)-apply(X2,2,median),ylab='Mean-Median')
    abline(h=1)
    #Call Identify 
    if(ident==T){
      identify(apply(X2,2,mean)-apply(X2,2,median),labels=paste(1:length(X[1,]),colnames(X)))
    }
  }
  
  #Plot sd
  if(sum(type==2)>0){
    plot(apply(X,2,sd),ylab='Standard Dev.')
    #Call Identify 
    if(ident==T){
      identify(apply(X,2,sd),labels=paste(1:length(X[1,]),colnames(X)))
    }
  }
  
  #Plot skewness
  if(sum(type==3)>0){
    if(corr==T){par(mfrow=c(2,1))}
    plot(apply(X,2,skewness),ylab='Skewness')
    if(SES==T){abline(h=c(momentsoutlier(X)$SES*3,-momentsoutlier(X)$SES*3))}
    #Call Identify 
    if(ident==T){
      identify(apply(X,2,skewness),labels=paste(1:length(X[1,]),colnames(X)))
    }
    #Plot X vs log(X)
    if(corr==T){
      X2=X
      X2[,momentsoutlier(X)$S.ind]=log(X2[,momentsoutlier(X)$S.ind]-min(X2[,momentsoutlier(X)$S.ind])+1)
      plot(apply(X2,2,skewness))
      abline(h=c(momentsoutlier(X)$SES*3,-momentsoutlier(X)$SES*3))
    }
  }
  
  #Plot kurtosis
  if(sum(type==4)>0){
    plot(apply(X,2,kurtosis),ylab='Kurtosis')
    if(SEK==T){abline(h=c(momentsoutlier(X)$SEK*3,-momentsoutlier(X)$SEK*3))}
    #Call Identify 
    if(ident==T){
      identify(apply(X,2,kurtosis),labels=paste(1:length(X[1,]),colnames(X)))
    }
  }
  #Restore par
  par(mfrow=c(1,1))
}

cplot<-function(X1,X2=NULL,setpar=NULL,ab=T,lml=F,xlab=NULL,ylab=NULL,...)
  #AKA combined plot, plots several columns against each other in 1 go
  #Input: matrices X1, X2; argument for par; plot abline?; plot lmline?; own x/y label
  #Output: plots columns of X1 against columns of X2
{
  #label of X and Y axis
  Xname=ifelse(is.null(xlab),deparse(substitute(X1)),xlab)
  Yname=ifelse(is.null(ylab),deparse(substitute(X2)),ylab)
  
  # also plot if only 1 matrix is specified
  if(is.null(X2)){X2 = matrix( 1:length(X1[,1]) ,nrow=length(X1[,1]),ncol=length(X1[1,]),byrow=F)}
  tmp = X1
  X1 = X2
  X2 = tmp
  
  #number of X1-variables
  d=length(X1[1,])
  
  #Determine plot par, if d>4 you should choose own par 
  # unless sqrt(d) is integer
  if(d<4){par(mfrow=c(d,1))}
  if(round(sqrt(d))==sqrt(d)){par(mfrow=c(sqrt(d),sqrt(d)))}
  if(!is.null(setpar)){par(mfrow=setpar)}
  if(!(d<4)&&!(round(sqrt(d))==sqrt(d))&&is.null(setpar)){stop('set a par')}
  
  #Plot pairwise columns of X1 and X2, possibly with lines
  for(i in 1:d){
    plot(X1[,i],X2[,i],xlab=paste(Xname,i),ylab=paste(Yname,i),...)
    Cf=lm(X2[,i]~X1[,i])$coef
    if(lml){
      abline(Cf,col=2)}
    if(ab){
      abline(0,sign(Cf[2]),col=3,lty=3)
    }
  }
  par(mfrow=c(1,1))
}


pplot<-function(X,Y,ab=T,setpar=NULL,Yvar=F,...)
  #AKA prediction plot, plots columns of X against Y (one should be vector)
  #Input: matrices X, Y; argument for par; plot abline?; own x/y label, multiple Y-variables?
  #Output: plots columns of X against Y or vice versa 
{
  # dimension of X or Y
  d=ifelse(Yvar,length(Y[1,]),length(X[1,]))
  
  #set a par
  if(d<4){par(mfrow=c(d,1))}
  if(round(sqrt(d))==sqrt(d)){par(mfrow=c(sqrt(d),sqrt(d)))}
  if(!(d<4)&&!(round(sqrt(d))==sqrt(d))&&!is.null(setpar)){par(mfrow=setpar)}
  if(!(d<4)&&!(round(sqrt(d))==sqrt(d))&&is.null(setpar)){stop('set a par')}
  
  #plot X[,i] vs Y or X vs Y[,i]
  for(i in 1:d){
    if(!Yvar){plot(X[,i],Y,...);if(ab){abline(lm(Y~X[,i])$coef,col=2)}}
    if(Yvar){plot(X,Y[,i],...);if(ab){abline(lm(Y[,i]~X[,1])$coef,col=2)}}
  }
  par(mfrow=c(1,1))
}

error.bar <- function(x, y, upper, lower=-upper,signif=T,...)
  #input: x (usually 1:p), y (usually loadings), upper and lower CI vectors, plot only significant y?
  #output:  plot of CI intervals for each (signif) y
{
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  incl=x
  if(signif){incl=which(lower*upper>=0)}
  x2=x[incl]
  y2=y[incl]
  lower2=lower[incl]
  upper2=upper[incl]
  plot(x,y,ylim=c(min(c(lower,upper)),max(c(lower,upper))),type='n',...)
  arrows(x2,upper2, x2, lower2, angle=90, code=3, length=0.1)
  abline(h=0,col=2,lwd=1,lty=1)
}

findr <- function(lambda,Cxt.=bla$a,Cxto.=bla$b,Ctt.=bla$c,Ctto.=bla$d,Ctoto.=bla$e,type=1)
{
  Cxt = Cxt.
  Cxto = Cxto.
  Ctt = c(Ctt.)
  Ctto = c(Ctto.)
  Ctoto = c(Ctoto.)
  
  alp1 = Ctt+lambda[1]
  alp2 = Ctoto+lambda[2]
  
  Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto^2)
  Wl = (Cxt - Pl*Ctto) / alp1
  
  tst1=crossprod(Wl)
  tst2=crossprod(Pl)
  if(type==1){return(c(tst1-1,tst2-1))}
  return(list(P=Pl,W=Wl))
}

findr_uni <- function(lambda2,Cxt.=bla$a,Cxto.=bla$b,Ctt.=bla$c,Ctto.=bla$d,Ctoto.=bla$e,first_root=T,type=1)
{
  Cxt = Cxt.
  Cxto = Cxto.
  Ctt = c(Ctt.)
  Ctto = c(Ctto.)
  Ctoto = c(Ctoto.)
  
  alp2 = Ctoto+lambda2
  
  a_rt = ssq(Cxto) - alp2^2
  b_rt = 2*Ctt*ssq(Cxto) - 2*(crossprod(Cxto,Cxt)*Ctto + Ctt*alp2^2 - alp2*Ctto^2)
  c_rt = ssq(Cxto)*Ctt^2 - 2*crossprod(Cxto,Cxt)*Ctt*Ctto + ssq(Cxt)*Ctto^2 - 
    Ctt^2*alp2^2 + 2*Ctt*alp2*Ctto^2 - Ctto^4
  
  D_rt = b_rt^2 - 4*a_rt*c_rt
  #if(D_rt < 0){stop("Discriminant negative!")}
  lambda1 = (-b_rt - (2*first_root-1)*sqrt(D_rt))/(2*a_rt)
  alp1 = c(Ctt+lambda1)
  
  Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto^2)
  Wl = (Cxt - Pl*Ctto) / alp1
  
  tst1=crossprod(Wl)
  tst2=crossprod(Pl)
  if(type==1){return(tst1-1)}
  return(list(P=Pl,W=Wl,l=lambda1))
}

jacr = function(lambda,Cxt.=bla$a,Cxto.=bla$b,Ctt.=bla$c,Ctto.=bla$d,Ctoto.=bla$e,type=1)
{
  Cxt = Cxt.
  Cxto = Cxto.
  Ctt = c(Ctt.)
  Ctto = c(Ctto.)
  Ctoto = c(Ctoto.)
  
  alp1 = Ctt+lambda[1]
  alp2 = Ctoto+lambda[2]
  Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto^2)
  Wl = (Cxt - Pl*Ctto) / alp1
  
  Pteller = (Cxto*alp1 - Cxt*Ctto)
  Pnoemer = (alp1*alp2 - Ctto^2)
  
  P1 = (Pnoemer * Cxto - Pteller * alp2) / Pnoemer^2
  P2 = -Pteller * alp1 / Pnoemer^2
  W1 = (alp1 * -Ctto*P1 - (Cxt-Pl*Ctto)) / alp1^2
  W2 = -Ctto/alp1*P2
  
  return(matrix(c(2*t(Wl)%*%W1,2*t(Pl)%*%P1,2*t(Wl)%*%W2,2*t(Pl)%*%P2),nrow=2))  
}

EMl <- function(W.=W, C.=C, P_Yosc.=PYo, P_Xosc.=PXo, B_T.=B_T, X.=X, Y.=Y,sigX.=sigX,sigY.=sigY,sigH.=sigH,sigT.=sigT,sigTo.=sigTo,sigUo.=sigUo,sr1=0*1:2,sr2=0*1:2,debug=FALSE,debug2=FALSE)
{
  N = length(X.[,1])
  p = length(X.[1,])
  q = length(Y.[1,])
  a = length(W.[1,])
  nx= length(P_Yosc.[1,])
  ny= length(P_Xosc.[1,])
  
  SX. = sigT.^2*tcrossprod(W.)+sigTo.^2*tcrossprod(P_Yosc.)+sigX.^2*diag(1,nrow=p)
  SY. = sigT.^2*tcrossprod(C.%*%t(B_T.))+tcrossprod(C.)*sigH.^2+sigUo.^2*tcrossprod(P_Xosc.)+sigY.^2*diag(1,nrow=q)
  SXY. = sigT.^2*(W.%*%B_T.%*%t(C.))
  #invSX. = woodinv(sigX.,cbind(W.,P_Yosc.))
  #invSY. = woodinv(sigY.,cbind(C.,P_Xosc.))
  Sfull = blockm(SX.,SXY.,SY.)
  invSfull = solve(Sfull)
  invSX. = solve(SX.)
  invSY. = solve(SY.)
  
  Cxx = 1/N*crossprod(X.)
  Cxy = 1/N*crossprod(X.,Y.)
  Cyy = 1/N*crossprod(Y.)
  Cxxyy = 1/N*crossprod(cbind(X.,Y.))
  
  Cxt = sigT.^2* (cbind(Cxx,Cxy)%*%invSfull%*%rbind(W.,C.%*%B_T.))
  Cxto = sigTo.^2* (Cxx%*%invSX.%*%P_Yosc.)
  Ctto = sigT.^2*sigTo.^2* (-t(rbind(W.,C.%*%B_T.))%*%invSfull%*%rbind(P_Yosc.,matrix(0,q)) + 
                              t(rbind(W.,C.%*%B_T.))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(P_Yosc.,matrix(0,q)))
  Ctt = sigT.^2* (diag(1,a) - sigT.^2*t(rbind(W.,C.%*%B_T.))%*%invSfull%*%rbind(W.,C.%*%B_T.) + 
                    sigT.^2*t(rbind(W.,C.%*%B_T.))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(W.,C.%*%B_T.))
  Ctoto = sigTo.^2* (diag(1,nx) - sigTo.^2*t(rbind(P_Yosc.,matrix(0,q)))%*%invSfull%*%rbind(P_Yosc.,matrix(0,q)) + 
                       sigTo.^2*t(rbind(P_Yosc.,matrix(0,q)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(P_Yosc.,matrix(0,q)))
  
  Cyt = sigT.^2* (cbind(t(Cxy),Cyy)%*%invSfull%*%rbind(W.,C.%*%B_T.))
  
  if(debug==T){return(list(a=Cxt,b=Cxto,c=Ctt,d=Ctto,e=Ctoto))}
  
  #Wnew = (Cxt - Cxto%*%solve(Ctoto)%*%t(Ctto))%*%solve(Ctt - Ctto%*%solve(Ctoto)%*%t(Ctto))
  #PYonew = (Cxto - Cxt%*%solve(Ctt)%*%(Ctto))%*%solve(Ctoto - t(Ctto)%*%solve(Ctt)%*%(Ctto))
  
  #############################
  #fr=list(root = c(0,0))
  #sr1 = c(1/c(crossprod(W.))*(t(W.)%*%Cxt - t(W.)%*%W.%*%Ctt - t(W.)%*%P_Yosc.%*%Ctto),
  #        1/c(crossprod(P_Yosc.))*(t(P_Yosc.)%*%Cxto - t(P_Yosc.)%*%W.%*%Ctto - t(P_Yosc.)%*%P_Yosc.%*%Ctoto))
  fr = multiroot(findr,start=sr1,maxit=1e6,jacfunc=jacr,Cxt.=Cxt,Cxto.=Cxto,Ctt.=Ctt,Ctto.=Ctto,Ctoto.=Ctoto)
  
  Wnew = findr(fr$root,Cxt.=Cxt,Cxto.=Cxto,Ctt.=Ctt,Ctto.=Ctto,Ctoto.=Ctoto,type=2)$W
  PYonew = findr(fr$root,Cxt.=Cxt,Cxto.=Cxto,Ctt.=Ctt,Ctto.=Ctto,Ctoto.=Ctoto,type=2)$P
  if(abs(vnorm(Wnew)+vnorm(PYonew)-2)>1e-6){
    if(abs(vnorm(Wnew)-1)>1e-6){Wnew = orth(findr(c(0,0),Cxt.=Cxt,Cxto.=Cxto,Ctt.=Ctt,Ctto.=Ctto,Ctoto.=Ctoto,type=2)$W)}
    if(abs(vnorm(PYonew)-1)>1e-6){PYonew = orth(findr(c(0,0),Cxt.=Cxt,Cxto.=Cxto,Ctt.=Ctt,Ctto.=Ctto,Ctoto.=Ctoto,type=2)$P)}
    warning('X:forced unit norm')
  }
  #############################
  
  Cyu = sigT.^2* (cbind(t(Cxy),Cyy)%*%invSfull%*%rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))
  Cyuo = sigTo.^2* (Cyy%*%invSY.%*%P_Xosc.)
  Cuuo = sigT.^2*sigUo.^2* (-t(rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))%*%invSfull%*%rbind(matrix(0,p),P_Xosc.) + 
                              t(rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(matrix(0,p),P_Xosc.))
  Cuu = sigT.^2* (crossprod(B_T.) - sigT.^2*t(rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))%*%invSfull%*%rbind(W.%*%B_T.,C.%*%crossprod(B_T.)) + 
                    sigT.^2*t(rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(W.%*%B_T.,C.%*%crossprod(B_T.)))
  Cuouo = sigUo.^2*(diag(1,ny) - sigUo.^2*t(rbind(matrix(0,p),P_Xosc.))%*%invSfull%*%rbind(matrix(0,p),P_Xosc.) + 
                      sigUo.^2*t(rbind(matrix(0,p),P_Xosc.))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(matrix(0,p),P_Xosc.))
  
  #Cnew = (Cyu - Cyuo%*%solve(Cuouo)%*%t(Cuuo))%*%solve(Cuu - Cuuo%*%solve(Cuouo)%*%t(Cuuo))
  #PXonew = (Cyuo - Cyu%*%solve(Cuu)%*%(Cuuo))%*%solve(Cuouo - t(Cuuo)%*%solve(Cuu)%*%(Cuuo))
  
  if(debug2==T){return(list(a=Cyu,b=Cyuo,c=Cuu,d=Cuuo,e=Cuouo))}
  
  #############################
  #fr2 = list(root = c(0,0))
  #sr2 = c(1/c(crossprod(C.))*(t(C.)%*%Cyu - t(C.)%*%C.%*%Cuu - t(C.)%*%P_Xosc.%*%Cuuo),
  #        1/c(crossprod(P_Xosc.))*(t(P_Xosc.)%*%Cyuo - t(P_Xosc.)%*%C.%*%Cuuo - t(P_Xosc.)%*%P_Xosc.%*%Cuouo))
  fr2 = multiroot(findr,start=sr2,maxit=1e6,jacfunc=jacr,Cxt.=Cyu,Cxto.=Cyuo,Ctt.=Cuu,Ctto.=Cuuo,Ctoto.=Cuouo)
  
  Cnew = findr(fr2$root,Cxt.=Cyu,Cxto.=Cyuo,Ctt.=Cuu,Ctto.=Cuuo,Ctoto.=Cuouo,type=2)$W
  PXonew = findr(fr2$root,Cxt.=Cyu,Cxto.=Cyuo,Ctt.=Cuu,Ctto.=Cuuo,Ctoto.=Cuouo,type=2)$P
  if(abs(vnorm(Cnew)+vnorm(PXonew)-2)>1e-6){
    if(abs(vnorm(Cnew)-1)>1e-6){Cnew = orth(findr(c(0,0),Cxt.=Cyu,Cxto.=Cyuo,Ctt.=Cuu,Ctto.=Cuuo,Ctoto.=Cuouo,type=2)$W)}
    if(abs(vnorm(PXonew)-1)>1e-6){PXonew = orth(findr(c(0,0),Cxt.=Cyu,Cxto.=Cyuo,Ctt.=Cuu,Ctto.=Cuuo,Ctoto.=Cuouo,type=2)$P)}
    warning('Y:forced unit norm')
  }
  #############################
  
  Bnew = vnorm(c(Cyt))/vnorm(c(Cxt))
  
  Cee = diag(sigX.^2,p) - t(rbind(diag(sigX.^2,p),diag(0,q,p)))%*%invSfull%*%rbind(diag(sigX.^2,p),diag(0,q,p)) + 
    t(rbind(diag(sigX.^2,p),diag(0,q,p)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(diag(sigX.^2,p),diag(0,q,p))
  sigXnew = sqrt(c(sum(diag(Cee))/p))
  
  Cff = diag(sigY.^2,q) - t(rbind(diag(0,p,q),diag(sigY.^2,q)))%*%invSfull%*%rbind(diag(0,p,q),diag(sigY.^2,q)) + 
    t(rbind(diag(0,p,q),diag(sigY.^2,q)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(diag(0,p,q),diag(sigY.^2,q))
  sigYnew = sqrt(c(sum(diag(Cff))/q))
  
  Chh = diag(sigH.^2,1) - t(rbind(matrix(0,p),sigH.^2*C.))%*%invSfull%*%rbind(matrix(0,p),sigH.^2*C.) + 
    t(rbind(matrix(0,p),sigH.^2*C.))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(matrix(0,p),sigH.^2*C.)
  sigHnew = sqrt(c(Chh))
  
  siglatnew = sqrt(c(Ctt,Ctoto,Cuouo))
  
  #print(c(Cut%*%solve(Ctt)))
  #list(What=Wnew/vnorm(Wnew),PYohat=PYonew/vnorm(PYonew),Chat=Cnew/vnorm(Cnew),PXohat=PXonew/vnorm(PXonew),Bhat=Bnew)
  list(What=Wnew,PYohat=PYonew,Chat=Cnew,PXohat=PXonew,Bhat=Bnew,sighat=c(sigXnew,sigYnew,sigHnew),siglathat=siglatnew,sr1.=sr1,sr2.=sr2)
}