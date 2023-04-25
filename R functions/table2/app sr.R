

app_sr=function(m=100,d=4,threshold=200,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){

  y0=mvrnorm(n=m,mu=mu0,Sigma=covar)
  y1=mvrnorm(n=m,mu=mu1,Sigma=covar)
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
    c=0
    r=0
    repeat{
        x=mvrnorm(n=1,mu=mu0,Sigma =covar)
        c=c+1
        print(c)
        r=(r+1)*dmvnorm(x,mean =estmu1,sigma =estcovar1)/dmvnorm(x,mean =estmu0,sigma =estcovar0)
        if(r>threshold){
          break
        }
        if(c>1000){
          break
        }
      }
  return(c)
}



app1_sr=function(m=100,d=4,threshold=200,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  c=0
  y0=mvrnorm(n=m,mu=mu0,Sigma=covar)
  y1=mvrnorm(n=m,mu=mu1,Sigma=covar)
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
    repeat{
      r=0
      c=0
      repeat{
        x=mvrnorm(n=1,mu=mu0,Sigma = covar)
        c=c+1
        print(c)
          r=(1+r)*dmvnorm(x,mean = estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
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
      if(false>=1){
        break
      }
    }
  if(false<1){
  nodetect=0
  repeat{
    x=mvrnorm(n=1,mu=mu1,Sigma = covar)
    c=c+1
    r=(1+r)*dmvnorm(x,mean = estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
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
  if(false>=1){
    return(c(NA,NA,false))
  }
}
