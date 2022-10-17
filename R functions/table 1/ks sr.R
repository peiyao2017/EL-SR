
library(ks)

ks_sr=function(m=100,threshold=1,mu0=0,mu1=0.5,covar=1){
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  wh0=hpi(y0,deriv.order =0)
  wh1=hpi(y1,deriv.order =0)
  f0=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  r=0
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
      r=(1+r)*f1(x)/f0(x)
      if(r>threshold){
        break
      }
      if(c>1000){
        break
      }
  }
  return(c)
}



ks1_sr=function(m=100,threshold=20,mu0=0,mu1=0.5,covar=1){
  false=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  wh0=hpi(y0,deriv.order =0)
  wh1=hpi(y1,deriv.order =0)
  f0=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    c=0
    r=0
      repeat{
        x=rnorm(n=1,mean=mu0,sd = covar)
        c=c+1
        print(c)
        r=(1+r)*f1(x)/f0(x)
        if(r>threshold){
          false=false+1
          break
        }
        print(c) 
        if(c>=49){
          break
        }
      }
    if(c>=49){
      break
    }
    if(false>=10){
      break
    }
  }
  if(false<10){
  nodetect=0
  repeat{
    x=rnorm(n=1,mean=mu1,sd = covar)
    c=c+1
    r=(r+1)*f1(x)/f0(x)
    if(r>threshold){
      break
    }
      if(c>=1000){
        nodetect=1
        break
      }
    }
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}

