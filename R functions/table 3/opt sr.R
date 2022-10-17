


opt_sr=function(d=2,threshold=200){
  c=0
  r=0
  repeat{
    x0=rgamma(d,shape = 9,rate=2)
    c=c+1
    print(c)
    r=(1+r)*prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
    if(r>threshold){
      break
    }
    if(c>1000){
      break
    }
  }
  return(c)
}


opt1_sr=function(d=2,threshold=200){
  false=0
  repeat{
    c=0
    r=0
    repeat{
      x0=rgamma(d,shape = 9,rate=2)
      c=c+1
      print(c)
      r=(1+r)*prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
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
      x0=rgamma(d,shape = 9,rate=2)+rep(0.5,times=d)
      c=c+1
      r=(1+r)*prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
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


