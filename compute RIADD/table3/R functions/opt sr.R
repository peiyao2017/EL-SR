


opt_sr=function(d=2,threshold=200,trunction=1000){
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
    if(c>trunction){
      break
    }
  }
  return(c)
}


opt1_sr=function(d=2,threshold=200){
  if(v<=1){
    repeat{
      c=0
      r=0
      repeat{
        x0=rgamma(d,shape = 9,rate=2)+mu1
        c=c+1
        print(c)
        r=(1+r)*prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
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
      x0=rgamma(d,shape = 9,rate=2)+rep(0.5,times=d)
      c=c+1
      r=(1+r)*prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
      if(r>threshold){
        break
      }
      if(c>1000){
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
    return(c(c-v,false,nodetect))
  }
  if(false>=1 ){
    return(c(NA,false,NA))
  }

  }
}

