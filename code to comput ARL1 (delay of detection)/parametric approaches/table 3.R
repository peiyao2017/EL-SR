


arl11=numeric()
arl12=numeric()
arl13=numeric()



repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    opt1_sr(threshold =7)
  }
  arl11=c(arl11,a[a!="false"])
  if(length(arl11)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    app1_sr(threshold =3005.5,m=20)
  }
  arl12=c(arl12,a[a!="false"])
  if(length(arl12)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    library(ks)
    ks1_sr(threshold =192246.5,m=20)
  }
  arl13=c(arl13,a[a!="false"])
  if(length(arl13)>2000){
    break
  }
}






x=data.frame(arl11=arl11[1:2000],arl12=arl12[1:2000],arl13=arl13[1:2000])

for(i in 1:3){
  x[,i]=as.numeric(x[,i])
}

delay=matrix(0,nrow=3,ncol=1000)
nodetect=matrix(0,nrow=3,ncol=1000)
for(i in 1:3){
  for(j in 1:1000){
    delay[i,j]=x[2*j-1,i]
    nodetect[i,j]=x[2*j,i]
  }
}


mean(delay[1,])
sd(delay[1,])/sqrt(1000)
sum(nodetect[1,])


mean(delay[2,])
sd(delay[2,])/sqrt(1000)
sum(nodetect[2,])


mean(delay[3,])
sd(delay[3,])/sqrt(1000)
sum(nodetect[3,])






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


ks1_sr=function(m=20,d=2,threshold=200,mu1=rep(0.5,times=d)){
  false=0
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  wh0=numeric()
  wh1=numeric()
  for(i in 1:d){
    wh0[i]=hpi(y0[,i],deriv.order =0)
    wh1[i]=hpi(y1[,i],deriv.order =0)
  }
  f0=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    r=0
    c=0
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)
      c=c+1
      print(c)
      r=(1+r)*f1(x)/f0(x)
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
      r=(r+1)*f1(x)/f0(x)
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
