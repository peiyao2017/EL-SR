




app_sr=function(m=100,d=4,threshold=200,lambda0=4,lambda1=4.5){
  c=0
  r=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rpois(n=d,lambda=lambda0)
    y1[i,]=rpois(n=d,lambda=lambda1)
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    x=rpois(n=d,lambda=lambda0)
    c=c+1
    print(c)
    r=(1+r)*dmvnorm(x,mean=estmu1,sigma = estcovar1)/dmvnorm(x,mean=estmu0,sigma = estcovar0)
        
      if(r>threshold){
        break
      }
    if(c>1000){
      break
    }
  }
  return(c)
}



app1_sr=function(m=100,d=4,threshold=200,lambda0=4,lambda1=4.5){
  false=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rpois(n=d,lambda=lambda0)
    y1[i,]=rpois(n=d,lambda=lambda1)
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    c=0
    r=0
    repeat{
      x=rpois(n=d,lambda=lambda0)
      c=c+1
      print(c)
      r=(1+r)*dmvnorm(x ,mean=estmu1,sigma = estcovar1)/dmvnorm(x ,mean=estmu0,sigma = estcovar0)
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
    x=rpois(n=d,lambda=lambda1)
    c=c+1
    r=(1+r)*dmvnorm(x ,mean=estmu1,sigma = estcovar1)/dmvnorm(x ,mean=estmu0,sigma = estcovar0)
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
