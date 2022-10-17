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





app1_sr=function(m=100,threshold=50,mu0=0,mu1=0.5,covar=1){
  false=0
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    c=0
    r=0
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
      r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
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
  nodetect=0
  if(false<10){
    repeat{
      x=rnorm(n=1,mean=mu1,sd = covar)
      c=c+1
      r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
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




ks1_sr=function(m=100,threshold=20,mu0=0,mu1=0.5,covar=1){
  false=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  wh0=hpi(y0,deriv.order =0)
  wh1=hpi(y1,deriv.order =0)
  f0=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    c=0
    r=0
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
      r=(1+r)*f1(x)/f0(x)
      if(r>threshold){
        false=false+1
        break
      }
      print(c) 
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
      x=rnorm(n=1,mean=mu1,sd = covar)
      c=c+1
      r=(r+1)*f1(x)/f0(x)
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


opt1_sr=function(threshold=1,mu0=0,mu1=0.5,covar=1){
  false=0
  repeat{
    c=0
    r=0
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
      r=(1+r)*dnorm(x,mean=mu1,sd=covar)/dnorm(x,mean=mu0,sd=covar)
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
      x=rnorm(n=1,mean=mu1,sd = covar)
      c=c+1
      r=(1+r)*dnorm(x,mean=mu1,sd=covar)/dnorm(x,mean=mu0,sd=covar)
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






