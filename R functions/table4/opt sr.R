opt_sr=function(d=4,threshold=200){
  c=0
  r=0
  repeat{
      x0=rgamma(d,shape = 10,rate=2)
      c=c+1
      r=(1+r)*prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
      if(r>threshold){
        break
      }
      if(c>1000){
        break
      }
  }
  return(c)
}




opt1_sr=function(d=4,threshold=200){
  false=0
  repeat{
    c=0
    r=0
    repeat{
        x0=rgamma(d,shape = 10,rate=2)
        c=c+1
        print(c)
        r=(1+r)*prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
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
    x0=rgamma(d,shape = 5,rate=1)
    c=c+1
    print(c)
    r=(r+1)*prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
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
  if(false>=1 ){
    return(c(NA,NA,false))
  }
}