
EMstep_W <- function(W.=W, C.=C, B_T.=B_T, X.=X, Y.=Y,
                     sigX.=sigX,sigY.=sigY,sigH.=sigH,sigT.=sigT,
                     debug=FALSE,debug2=FALSE,invS = t(0))
{
  N = nrow(X.)
  p = ncol(X.)
  q = ncol(Y.)
  a = dim(W.)[2]
  B_T.=c(B_T.)
  sigT. = c(sigT.)
  
  XY = cbind(X.,Y.)
#   SX. = sigT.^2*tcrossprod(W.)+sigX.^2*diag(1,nrow=p)
#   SY. = sigT.^2*B_T.^2*tcrossprod(C.)+tcrossprod(C.)*sigH.^2+sigY.^2*diag(1,nrow=q)
#   SXY. = sigT.^2*B_T.*(W.%*%t(C.))
#   #  invSX. = woodinv(sigX.,cbind(W.))
#   #  invSY. = woodinv(sigY.,cbind(C.))
#   Sfull = blockm(SX.,SXY.,SY.)
#   invSfull = solve(Sfull)
#   
#   
#   Cxx = crossprod(X.)/N
#   Cxy = crossprod(X.,Y.)/N
#   Cyy = crossprod(Y.)/N
#   Cxxyy = crossprod(XY)/N
#   
#   covXYT = rbind(sigT.^2*W. , sigT.^2*C.*B_T.)
#   del_T = invSfull %*% covXYT
#   Del_T = sigT.^2 - t(covXYT)%*%invSfull%*%covXYT
#   
#   Cxt = cbind(Cxx,Cxy)%*%del_T
#   Ctt = t(del_T)%*%Cxxyy%*%del_T + Del_T
#   #  Cyt = sigT.^2* (cbind(t(Cxy),Cyy)%*%invSfull%*%rbind(W.,C.%*%B_T.))
#   
#   if(debug==T){return(list(a=Cxt,b=0*Cxt,c=Ctt,d=0,e=0))}
#   
#   covXYU = rbind(sigT.^2*B_T.*W. , (sigT.^2*B_T.^2+sigH.^2)*C.)
#   del_U = invSfull %*% covXYU
#   Del_U = sigT.^2*B_T.^2+sigH.^2 - t(covXYU)%*%invSfull%*%covXYU
#   
#   Cyu = cbind(t(Cxy),Cyy)%*%del_U
#   Cuu = t(del_U)%*%Cxxyy%*%del_U + Del_U
#   
#   if(debug2==T){return(list(a=Cyu,b=0*Cyu,c=Cuu,d=0,e=0))}
#   
#   Cut = sigT.^2*B_T. - t(covXYU)%*%invSfull%*%covXYT + t(del_U)%*%Cxxyy%*%del_T
#   Bnew = Cut / c(Ctt)
#   
#   Cee = diag(sigX.^2,p) - t(rbind(diag(sigX.^2,p),diag(0,q,p)))%*%invSfull%*%rbind(diag(sigX.^2,p),diag(0,q,p)) + 
#     t(rbind(diag(sigX.^2,p),diag(0,q,p)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(diag(sigX.^2,p),diag(0,q,p))
#   sigXnew = sqrt(tr(Cee)/p)
#   
#   Cff = diag(sigY.^2,q) - t(rbind(diag(0,p,q),diag(sigY.^2,q)))%*%invSfull%*%rbind(diag(0,p,q),diag(sigY.^2,q)) + 
#     t(rbind(diag(0,p,q),diag(sigY.^2,q)))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(diag(0,p,q),diag(sigY.^2,q))
#   sigYnew = sqrt(tr(Cff)/q)
#   
#   Chh = diag(sigH.^2,1) - t(rbind(matrix(0,p),sigH.^2*C.))%*%invSfull%*%rbind(matrix(0,p),sigH.^2*C.) + 
#     t(rbind(matrix(0,p),sigH.^2*C.))%*%invSfull%*%Cxxyy%*%invSfull%*%rbind(matrix(0,p),sigH.^2*C.)
#   sigHnew = sqrt(c(Chh))
#   
  #Expec = list(Cxt,0*1:p,Ctt,0,0,Cyu,0*1:q,Cuu,0,0,Cut,diag(Cee),diag(Cff),Chh)
  Expec = EMstepC(W.,C.,0*W.,0*C.,t(B_T.),XY,
                  sigX.^2,sigY.^2,t(sigH.)^2,t(sigT.)^2,t(0),t(0),invS)
  #return(list(C = Expec, R = list(Cxt,0*1:p,Ctt,0,0,Cyu,0*1:q,Cuu,0,0,Cut,diag(Cee),diag(Cff),Chh)))
  Cxt = Expec$Cxt
  #Cxto = Expec$Cxto
  Ctt = Expec$Ctt
  #Ctto = Expec$Ctto
  #Ctoto = Expec$Ctoto
  
  Cyu = Expec$Cyu
  #Cyuo = Expec$Cyuo
  Cuu = Expec$Cuu
  #Cuuo = Expec$Cuuo
  #Cuouo = Expec$Cuouo
  
  Cut = Expec$Cut
  #############################
  
  
  Wnew = Cxt / vnorm(Cxt) #(Cxt - Cxto%*%solve(Ctoto)%*%t(Ctto))%*%solve(Ctt - Ctto%*%solve(Ctoto)%*%t(Ctto))
  #PYonew = P_Yosc. #(Cxto - Cxt%*%solve(Ctt)%*%(Ctto))%*%solve(Ctoto - t(Ctto)%*%solve(Ctt)%*%(Ctto))
  Cnew = Cyu / vnorm(Cyu) #(Cyu - Cyuo%*%solve(Cuouo)%*%t(Cuuo))%*%solve(Cuu - Cuuo%*%solve(Cuouo)%*%t(Cuuo))
  #PXonew = P_Xosc. #(Cyuo - Cyu%*%solve(Cuu)%*%(Cuuo))%*%solve(Cuouo - t(Cuuo)%*%solve(Cuu)%*%(Cuuo))
  Bnew = Cut / c(Ctt)
  sigXnew = sqrt(mean(Expec$diag_Cee))
  sigYnew = sqrt(mean(Expec$diag_Cff))
  sigHnew = sqrt(mean(Expec$diag_Cff))
  
  #   sigXnew = sqrt(c( tr(Cxx - 2*Cxt%*%t(Wnew) + Wnew%*%Ctt%*%t(Wnew)) )/p)
  #   sigYnew = sqrt(c( tr(Cyy - 2*Cyu%*%t(Cnew) + Cnew%*%Cuu%*%t(Cnew)) )/q)
  #  sigHnew = 1#sqrt(c(Cuu - 2*Cut*B_T. + Ctt*B_T.^2))
  siglatnew = c(c(sigHnew,sqrt(Ctt)))
  #############################
  
  list(What=Wnew,Chat=Cnew,Bhat=matrix(Bnew),sighat=c(sigXnew,sigYnew),siglathat=siglatnew)
  #  list(What=Wnew,Chat=Cnew,Bhat=B_T.,sighat=c(sigX.,sigY.,sigH.),siglathat=sigT.)
}

PPLSi <- function(X,Y,EMsteps=20,atol=1e-4,ret=c("o2m","random","uniform"),critfunc=function(x){x},print_logvalue=FALSE)
{
  suppressWarnings({if(ret == "o2m"){
    sim=o2m(X,Y,1,0,0)
    Wnw=sim$W.;Cnw=sim$C.;Bnw=sim$B_T.;
    signw=sqrt(c(sum(diag(crossprod(sim$E)))/N/p,sum(diag(crossprod(sim$F)))/N/q));
    siglatnw=sqrt(c(ssq(sim$H_U),ssq(sim$Tt))/N)
    #set.seed(1);Wnw=orth(runif(p));PYonw=orth(runif(p));Cnw=orth(runif(q));PXonw=orth(runif(q)); Bnw = matrix(rchisq(1,1)); siglatnw =  rchisq(4,100)/100; signw = rchisq(2,10)/100
    #Wnw=orth(1:p);PYonw=orth(1:p);Cnw=orth(1:q);PXonw=orth(1:q);
    ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }})
  suppressWarnings({if(ret == "random"){
    Wnw=orth(runif(p));Cnw=orth(runif(q)); Bnw = matrix(rchisq(1,1)); siglatnw =  rchisq(2,100)/100; signw = rchisq(2,10)/100
    #Wnw=orth(1:p);PYonw=orth(1:p);Cnw=orth(1:q);PXonw=orth(1:q);
    ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }})
  suppressWarnings({if(ret == "uniform"){
    Wnw=orth(rep(1,p));Cnw=orth(rep(1,q)); Bnw = matrix(1); siglatnw =  c(1,1); signw = c(1/p,1/q)
    ret = list(W=Wnw,C=Cnw,B=Bnw,sig=c(signw,siglatnw))
  }})
  p = length(ret$W); q = length(ret$C); 
  Dat = cbind(X,Y)
  Wnw=ret$W;Cnw=ret$C;Bnw=c(ret$B);signw=ret$sig[1:2];siglatnw=ret$sig[3:4]
  
  logvalue=1:(EMsteps+1)*NA
  
  loglClist=loglC(Wnw,Cnw,0*Wnw,0*Cnw,t(Bnw),Dat,signw[1]^2,signw[2]^2,t(siglatnw[1]^2),t(siglatnw[2]^2),t(0),t(0))
  logvalue[1] = loglClist$val
  
  for(i in 1:EMsteps){
    fit=EMstep_W(Wnw,Cnw,Bnw,X,Y,
                 signw[1],signw[2],siglatnw[1],siglatnw[2],0,0,loglClist$invS)
    Bnw = (fit$Bhat)
    Wnw = (fit$W)
    PYonw = (fit$PY)
    Cnw = (fit$C)
    PXonw = (fit$PX)
    signw = fit$sigh
    siglatnw = fit$sigl
    loglClist=loglC(Wnw,Cnw,0*Wnw,0*Cnw,t(Bnw),Dat,signw[1]^2,signw[2]^2,t(siglatnw[1]^2),t(siglatnw[2]^2),t(0),t(0))
    logvalue[i+1] = loglClist$val
    if(critfunc(logvalue[i+1]-logvalue[i])<atol){break}
    on.exit(return(list( W=c(Wnw),C=c(Cnw),B=Bnw,sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = logvalue[i]-logvalue[i-1], Number_steps = i)))
  }
  last_incr = logvalue[i+1]-logvalue[i]
  if(any(diff(logvalue[0:i+1])<0)){warning("Not monotone")}
  if(critfunc(logvalue[i+1]-logvalue[i])>=atol){warning(paste("Not converged, last increment was",critfunc(logvalue[i+1]-logvalue[i])))}
  if(!print_logvalue){logvalue=NA}
  ret2 = list( W=c(Wnw),C=c(Cnw),B=Bnw,sig=c(signw,siglatnw),logvalue=logvalue[0:i+1],Last_increment = last_incr, Number_steps = i)
  return(ret2)
}

PPLS <- function(X,Y,nr_comp=a,EMsteps=20,atol=1e-4,ret=c("o2m","random","uniform"),critfunc=function(x){x},print_logvalue=FALSE)
{
  a = nr_comp
  Wn = matrix(NA,p,a)
  Cn = matrix(NA,q,a)
  Bn = sigXn = sigYn = sigHn = sigTn = 1:a*NA
  other_output = list()
  Xc = X
  Yc = Y
  for(i in 1:a){
    fit = PPLSi(Xc,Yc,EMsteps,atol,ret,critfunc,print_logvalue=TRUE);
    Wn[,i] = c(fit$W);
    Cn[,i] = c(fit$C);
    Bn[i] = c(fit$B);
    sigXn[i] = fit$sig[1]
    sigYn[i] = fit$sig[2]
    sigHn[i] = fit$sig[3]
    sigTn[i] = fit$sig[4]
    Xc = Xc - Xc%*%tcrossprod(fit$W)
    Yc = Yc - Yc%*%tcrossprod(fit$C)
    other_output$Last_increment[i] = fit$Last
    other_output$Number_steps[i] = fit$Number
    other_output$Loglikelihoods[i] = logl_M(Wn[,1:i],Cn[,1:i],0*Wn,0*Cn,diag(Bn,i),cbind(X,Y),mean(sigXn[1:i]),mean(sigYn[1:i]),mean(sigHn[1:i]),diag(sigTn[1:i],i),diag(0,a),diag(0,a))
  }
  ret2 = list( W=Wn,C=Cn,B=Bn,sig=cbind(sigX=sigXn,sigY=sigYn,sigH=sigHn,sigT=sigTn),Other_output=other_output)
  
  return(ret2)
}
