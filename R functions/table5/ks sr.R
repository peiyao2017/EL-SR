




ks_sr=function(m=100,d=4,threshold=100,min0=-3,max0=3,min1=-3.5,max1=3.5){
  c=0
  r=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=runif(n=d,min=min0,max=max0)
    y1[i,]=runif(n=d,min=min1,max=max1)
    
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
    x=runif(n=d,min=min0,max=max0)
    c=c+1
    print(c)
    r=(r+1)*f1(x)/f0(x)
      if(r>threshold){
        break
      }
    if(c>1000){
      break
    }
  }
  return(c)
}



ks1_sr=function(m=100,d=4,threshold=100,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=runif(n=d,min=min0,max=max0)
    y1[i,]=runif(n=d,min=min1,max=max1)
    
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
      x=runif(n=d,min=min0,max=max0)
      c=c+1
      r=(r+1)*f1(x)/f0(x)
        if(r>threshold){
          false=false+1
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
  if(false<1 ){
  nodetect=0
  repeat{
    x=runif(n=d,min=min1,max=max1)
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
  return(c(c-50,nodetect,NA))
  }
  if(false>=1 ){
    return(c(NA,NA,false))
  }
}
