



app_sr=function(m=100,d=4,threshold=200,lambda0=4,lambda1=4.5,trunction=1000){
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
    if(c>trunction){
      break
    }
  }
  return(c)
}


app1_sr=function(m=100,d=4,threshold=200,lambda0=4,lambda1=4.5,v=0){
  if(v<=1){
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
        x=rpois(n=d,lambda=lambda1)
        c=c+1
        print(c)
        r=(1+r)*dmvnorm(x ,mean=estmu1,sigma = estcovar1)/dmvnorm(x ,mean=estmu0,sigma = estcovar0)
        if(r>threshold){
          break
        }
        if(c>=50){
          break
        }
      }
      if(c>=50){
        break
      }
      if(r>threshold ){
        break
      }
    }
    if(r>threshold){
      delay=c-v-1
      false=0
      nodetect=0
      return(c(delay,false,nodetect))
    }
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
      return(c(c-v-1,0,nodetect))
  }
  if(v>1){
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
  return(c(c-v-1,false,nodetect))
  }
  if(false>=1 ){
    return(c(NA,false,NA))
  }
  }
  
}


