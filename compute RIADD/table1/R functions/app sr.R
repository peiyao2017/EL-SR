

app_sr=function(m=100,threshold=100,mu0=0,mu1=0.5,covar=1,trunction=1000){
  c=0
  r=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
     r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
     print(r)
  if(r>threshold){
    break
  }
     if(c>trunction){
       break
     }
  }
  return(c)
}



app1_sr=function(m=100,threshold=50,mu0=0,mu1=0.5,covar=1,v=0){
  if(v<=1){
    c=0
    y0=rnorm(n=m,mean=mu0,sd=covar)
    y1=rnorm(n=m,mean=mu1,sd=covar)
    estmu0=mean(y0)
    estmu1=mean(y1)
    estcovar0=sd(y0)
    estcovar1=sd(y1)
    repeat{
      c=0
      r=0
      repeat{
        x=rnorm(n=1,mean=mu1,sd = covar)
        c=c+1
        print(c)
        r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
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
      if(r>threshold){
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
        x=rnorm(n=1,mean=mu1,sd = covar)
        c=c+1
        r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
        if(r>threshold){
          break
        }
        if(c>=1000){
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
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    c=0
    r=0
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
        r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
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
  nodetect=0
  if(false<1 ){
  repeat{
    x=rnorm(n=1,mean=mu1,sd = covar)
    c=c+1
    r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
    if(r>threshold){
      break
    }
    if(c>=1000){
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
