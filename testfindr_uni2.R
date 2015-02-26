source('C:/Users/selbouhaddani/OneDrive/LUMC/PhD/Rcode/Functions/findr_uni2.R')

bla=EMl(debug=T)

rootf = function(func,xstart,xrng,cl=NULL,maxiter=100,atol=1e-8,fact=10,gridnr=1000,...){
  tst1=tst2=12345
  xtst.min = xstart
  root1=root2=0
  
  K = min(maxiter,round(-log(atol)/log(fact))+1)
  
  for(i in 0:K){
    xtst = seq(xtst.min-xrng/(fact^i),xtst.min+xrng/(fact^i),length=gridnr)
    if(root2<K/10){tst1 = func(xtst,first_root=T,...)}
    if(root1<K/10){tst2 = func(xtst,first_root=F,...)}
    
    if(min(abs(tst1),na.rm=T)<=min(abs(tst2),na.rm=T)){
      xtst.min = xtst[order(abs(tst1))[1]]
      tst=tst1
      root1=root1+1
    }
    if(min(abs(tst2),na.rm=T)<min(abs(tst1),na.rm=T)){
      xtst.min = xtst[order(abs(tst2))[1]]
      tst=tst2
      root2=root2+1
    }
  }
  if(min(abs(tst),na.rm=T)>atol*10){warning(paste("min(abs(tst)) = ",min(abs(tst))))}
  list(fmin = func(xtst.min,first_root=(root1>root2),...),
       xmin = xtst.min,
       x_acc = xrng/(fact^i),
       y_acc = diff(range(tst)),
       first_root=(root1>root2))
}

