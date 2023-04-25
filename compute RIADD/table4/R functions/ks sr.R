




ks_sr=function(m=100,d=2,threshold=100,trunction=1000){
  c=0
  r=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 10,rate = 2)
    y1[i,]=rgamma(n=d,shape = 5,rate = 1)
    
  }
  wh0=numeric()
  wh1=numeric()
  for(i in 1:d){
    wh0[i]=hpi(y0[,i],deriv.order =0)
    wh1[i]=hpi(y1[,i],deriv.order =0)
  }
  f0=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    x=rgamma(n=d,shape = 10,rate = 2)
    c=c+1
    print(c)
    r=(r+1)*f1(x)/f0(x)
      if(r>threshold){
        break
      }
    if(c>trunction){
      break
    }
  }
  return(c)
}



ks1_sr=function(m=100,d=2,threshold=100,v=0){
 if(v<=1){
   y0=matrix(0,m,d)
   y1=matrix(0,m,d)
   for(i in 1:m){
     y0[i,]=rgamma(n=d,shape = 10,rate = 2)
     y1[i,]=rgamma(n=d,shape = 5,rate = 1)
     
   }
   wh0=numeric()
   wh1=numeric()
   for(i in 1:d){
     wh0[i]=hpi(y0[,i],deriv.order =0)
     wh1[i]=hpi(y1[,i],deriv.order =0)
   }
   f0=function(x=rep(0,times=d)){
     p=numeric()
     for(i in 1:m){
       p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
     }
     return(mean(p))
   }
   f1=function(x=rep(0,times=d)){
     p=numeric()
     for(i in 1:m){
       p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
     }
     return(mean(p))
   }
   repeat{
     c=0
     r=0
     repeat{
       x=rgamma(n=d,shape = 5,rate =1)
       c=c+1
       r=(r+1)*f1(x)/f0(x)
       if(r>threshold){
         break
       }
       if(c>=49){
         break
       }
     }
     if(c>=49){
       break
     }
     if(false>=1 ){
       break
     }
   }
   if(r>threshold){
     delay=c-v
     false=0
     nodetect=0
     return(c(delay,false,nodetect))
   }
     nodetect=0
     repeat{
       x=rgamma(n=d,shape = 5,rate = 1)
       c=c+1
       r=(1+r)*f1(x)/f0(x)
       if(r>threshold){
         break
       }
       if(c>1000){
         nodetect=1
         break
       }
     }
     if(c<=1000){
       delay=c-v
       false=0
       nodetect=0
       return(c(delay,false,nodetect))
     }
     if(c>1000){
       delay=c-v
       false=0
       nodetect=1
       return(c(delay,false,nodetect))
     }
   }
 if(v>1){
   false=0
   y0=matrix(0,m,d)
   y1=matrix(0,m,d)
   for(i in 1:m){
     y0[i,]=rgamma(n=d,shape = 10,rate = 2)
     y1[i,]=rgamma(n=d,shape = 5,rate = 1)
     
   }
   wh0=numeric()
   wh1=numeric()
   for(i in 1:d){
     wh0[i]=hpi(y0[,i],deriv.order =0)
     wh1[i]=hpi(y1[,i],deriv.order =0)
   }
   f0=function(x=rep(0,times=d)){
     p=numeric()
     for(i in 1:m){
       p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
     }
     return(mean(p))
   }
   f1=function(x=rep(0,times=d)){
     p=numeric()
     for(i in 1:m){
       p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
     }
     return(mean(p))
   }
   repeat{
     c=0
     r=0
     repeat{
       x=rgamma(n=d,shape = 10,rate = 2)
       c=c+1
       r=(r+1)*f1(x)/f0(x)
       if(r>threshold){
         false=false+1
         break
       }
       if(c>=v){
         break
       }
     }
     if(c>=v){
       break
     }
     if(false>=1 ){
       break
     }
   }
   if(false<1 ){
     nodetect=0
     repeat{
       x=rgamma(n=d,shape = 5,rate = 1)
       c=c+1
       r=(1+r)*f1(x)/f0(x)
       if(r>threshold){
         break
       }
       if(c>1000){
         nodetect=1
         break
       }
     }
     return(c(c-v,false,nodetect))
   }
   if(false>=1 ){
     return(c(NA,false,NA))
   }
 } 
}
