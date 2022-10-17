opt_sr=function(d=4,threshold=200,lambda0=4,lambda1=4.5){
  c=0
  r=0
  repeat{
      x0=rpois(n=d,lambda=lambda0)
      c=c+1
      r=(1+r)*prod(dpois(x0,lambda = lambda1)/dpois(x0,lambda = lambda0))
      if(r>threshold){
        break
      }
      if(c>1000){
        break
      }
  }
  return(c)
}




opt1_sr=function(d=4,threshold=200,lambda0=4,lambda1=4.5){
  false=0
  repeat{
    c=0
    r=0
    repeat{
        x0=rpois(n=d,lambda=lambda0)
        c=c+1
        print(c)
        r=(1+r)*prod(dpois(x0,lambda = lambda1)/dpois(x0,lambda = lambda0))
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
    x0=rpois(n=d,lambda=lambda1)
    c=c+1
    print(c)
    r=(r+1)*prod(dpois(x0,lambda = lambda1)/dpois(x0,lambda = lambda0))
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