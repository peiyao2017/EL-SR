
library(el.convex)



el_laplace=function(m=100,mu0=0,mu1=0.5,threshold=20){
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rep(0,times=m)
  x2=rep(0,times=m)
  y1=0.5
  y2=-0.5
  for(i in 1:m){
    x1[i]=exp(-y1*x[i])
    x2[i]=exp(-y2*x[i])
  }
  x3=mean(x1)
  x4=mean(x2)
  x5=c(x3,x4)
  x6=rnorm(n = m, mean=mu1,sd=1)
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1*x6[i]),exp(-y2*x6[i]))
  }
  x8=colMeans(x7)
  c=3
  x9=numeric()
  h0=matrix(0,c,2)
  h1=matrix(0,c,2)
  for(i in 1:c){
    x9[i]=rnorm(n=1,mean=mu0,sd=1)
  }
  for(i in 1:c){
    h0[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5
    h1[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8
  }
  z=numeric()
  repeat{
    for(i in 1:(c-2)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2)) 
      l1=el.test.newton(x=x11,mu=rep(0,2)) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
      z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
    }
    if(length(z)>0){
      print(sum(z))
      if(sum(z)<=threshold){
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
    }
    if(length(z)>0){
      if(sum(z)>threshold){
        break
      }
    }
    if(c>100){
      break
    }
  }
  if(length(z)==0){
    z=c(0)
  }
  if(sum(z)<=threshold){
    repeat{
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-2)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2)) 
        l1=el.test.newton(x=x11,mu=rep(0,2))
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
      if(length(z)>0){
        print(sum(z))
        if(sum(z)>threshold){
          print(sum(z))
          break
        }
      }
      if(length(z)>0){
        print(sum(z))
        if(sum(z)<=threshold){
          x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
          h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
        }
      }
      if(c>10000){
        break
      }
    }
  }
  return(c)
}



el1_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
  false=0
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rep(0,times=m)
  x2=rep(0,times=m)
  y1=0.5
  y2=-0.5
  for(i in 1:m){
    x1[i]=exp(-y1*x[i])
    x2[i]=exp(-y2*x[i])
  }
  x3=mean(x1)
  x4=mean(x2)
  x5=c(x3,x4)
  x6=rnorm(n = m, mean=mu1,sd=1)
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1*x6[i]),exp(-y2*x6[i]))
  }
  x8=colMeans(x7)
  repeat{
    c=3
    repeat{
      x9=numeric()
      h0=matrix(0,c,2)
      h1=matrix(0,c,2)
      for(i in 1:c){
        x9[i]=rnorm(n=1,mean=mu0,sd=1)
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5
        h1[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-2)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,2)) 
          l1=el.test.newton(x=x11,mu=rep(0,2)) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rnorm(n=1,mean=mu0,sd=1))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5)
          h1=rbind(h1,c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8)
        }
        if(length(z)>0){
          if(sum(z)>threshold){
            false=false+1
            break
          }
        }
        if(length(z)>0){
          if(sum(z)<=threshold){
            x9=rbind(x9,rnorm(n=1,mean=mu0,sd=1))
            c=c+1
            print(c)
            h0=rbind(h0,c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5)
            h1=rbind(h1,c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8)
          }
        }
        if(c>=49){
          break
        }
      }
      if(false>=1){
        break
      }
      if(c>=49){
        break
      }
    }
    if(c>=49){
      break
    }
    if(false>=1){
      break
    }
  }
  nodetect=0
  if(false<1){
  repeat{
    x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
    c=c+1
    h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
    h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
    z=numeric()
    for(i in 1:(c-2)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2)) 
      l1=el.test.newton(x=x11,mu=rep(0,2)) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    }
    if(length(z)>0){
      if(sum(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(sum(z)<=threshold){
        x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
    }
    if(c>=101){
      break
    }
  }
  if(length(z)==0){
    z=c(0)
  }
  if(sum(z)<=threshold){
    repeat{
      z=numeric()
      x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-2)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2)) 
        l1=el.test.newton(x=x11,mu=rep(0,2)) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
      if(length(z)>0){
        print(sum(z))
        if(sum(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        print(sum(z))
        if(sum(z)<=threshold){
          x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
          h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
        }
      }
      if(c>1000){
        nodetect=nodetect+1
        break
      }
    }
  }
  return(c(c-50,nodetect,false))
  }
  if(false){
    return(c(NA,NA,false))
  }
}



