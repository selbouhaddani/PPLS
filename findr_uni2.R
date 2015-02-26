findr_uni2 <- function(lambda2.,Cxt.=bla$a,Cxto.=bla$b,Ctt.=bla$c,Ctto.=bla$d,Ctoto.=bla$e,first_root=T,type=1)
{
  y=sapply(lambda2.,function(lambda2){
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
  
  D_rt = (b_rt^2 - 4*a_rt*c_rt) ; if(c(D_rt)<0){D_rt=NA}
  lambda1 = (-b_rt - (2*first_root-1)*sqrt(D_rt))/(2*a_rt)
  alp1 = c(Ctt+lambda1)
  
  Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto^2)
  Wl = (Cxt - Pl*Ctto) / alp1
  
  return(tst1=crossprod(Wl)-1)
  #tst2=crossprod(Pl)
  })
  return(y)
  #return((b_rt^2-2*b_rt+4*a_rt*c_rt)/(4*a_rt))
  #if(type==1){return(tst1-1)}
  #return(list(P=Pl,W=Wl,l=lambda1,xtop=-c(b_rt)/c(2*a_rt),ytop=c(b_rt^2-2*b_rt+4*a_rt*c_rt)/c(4*a_rt)))
}