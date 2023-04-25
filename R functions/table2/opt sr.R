

opt_sr=function(d=2,threshold=200,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  c=0
  r=0
  repeat{
    x=mvrnorm(n=1,mu=mu0,Sigma = covar)
    c=c+1
    print(c)
    r=(1+r)*dmvnorm(x,mean = mu1,sigma = covar)/dmvnorm(x,mean = mu0,sigma = covar)
      if(r>threshold){
        break
      }
    if(c>1000){
      break
    }
    }
  return(c)
}



opt1_sr=function(d=2,threshold=200,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  repeat{
    c=0
    r=0
    repeat{
    x=mvrnorm(n=1,mu=mu0,Sigma = covar)
    c=c+1
    print(c)
    r=(1+r)*dmvnorm(x,mean = mu1,sigma = covar)/dmvnorm(x,mean = mu0,sigma = covar)
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
    x=mvrnorm(n=1,mu=mu1,Sigma = covar)
    c=c+1
    r=(1+r)*dmvnorm(x,mean = mu1,sigma = covar)/dmvnorm(x,mean = mu0,sigma = covar)
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