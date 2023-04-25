opt_sr=function(d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  c=0
  r=0
  repeat{
      x0=runif(n=d,min=min0,max=max0)
      c=c+1
      r=(1+r)*prod(dunif(x0,min=min1,max=max1)/dunif(x0,min=min0,max=max0))
      if(r>threshold){
        break
      }
      if(c>1000){
        break
      }
  }
  return(c)
}




opt1_sr=function(d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  repeat{
    c=0
    r=0
    repeat{
        x0=runif(n=d,min=min0,max=max0)
        c=c+1
        print(c)
        r=(1+r)*prod(dunif(x0,min=min1,max=max1)/dunif(x0,min=min0,max=max0))
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
    x0=runif(n=d,min=min1,max=max1)
    c=c+1
    print(c)
    r=(r+1)*prod(dunif(x0,min=min1,max=max1)/dunif(x0,min=min0,max=max0))
    if(r>threshold){
      break
    }
      if(c>=1000){
        nodetect=1
        break
      }
  }
  return(c(c-50,nodetect,NA))
  }
  if(false>=10){
    return(c(NA,NA,false))
  }
}