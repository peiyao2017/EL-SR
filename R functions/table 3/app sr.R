




app_sr=function(m=20,d=2,threshold=200,mu1=rep(0.5,d)){
  c=0
  r=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    x=rgamma(n=d,shape = 9,rate = 2)
    c=c+1
    print(c)
    r=(r+1)*dmvnorm(x,mean = estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
    if(r>threshold){
      break
    }
    if(c>1000){
      break
    }
  }
  return(c)
}



app1_sr=function(m=20,d=2,threshold=200,mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    c=0
    r=0
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)
      c=c+1
      print(c)
      r=(r+1)*dmvnorm(x,mean =estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
      
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
    if(false>=10){
      break
    }
  }
  if(false<10){
    nodetect=0
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)+mu1
      c=c+1
      r=(r+1)*dmvnorm(x,mean =estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
      if(r>threshold){
        break
      }
      if(c>1000){
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


