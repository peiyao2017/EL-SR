


library(el.convex)



tael_laplace=function(m=50,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=200){
  x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(0.5,times=d)
  y2=rep(-0.5,times=d)
  for(i in 1:m){
    x1[i,]=exp(-y1%*%x[i,])
    x2[i,]=exp(-y2%*%x[i,])
  }
  x3=colMeans(x1)
  x4=colMeans(x2)
  x5=cbind(x3,x4)
  x6=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x6[i,]),exp(-y2%*%x6[i,]))
  }
  x8=colMeans(x7)
  c=d+1
  if(c<=100){
    x9=matrix(0,c,d)
    h0=matrix(0,c,2)
    h1=matrix(0,c,2)
    for(i in 1:c){
      x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
    }
    for(i in 1:c){
      h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
      h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
    }
    z=numeric()
    repeat{
      for(i in 1:(c-(d-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        if(i==c){
          x10s=x10
          x11s=x11
        }else{
          x10s=colMeans(x10)
          x11s=colMeans(x11)
        }
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        l0=l0$`-2LLR`
        l1=l1$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      if(c>100){
        break
      }
    }
  }
  if(sum(z)<=threshold){
    repeat{
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-(d-1))){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        if(i==nrow(current0)){
          x10s=x10
          x11s=x11
        }else{
          x10s=colMeans(x10)
          x11s=colMeans(x11)
        }
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        l0=l0$`-2LLR`
        l1=l1$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      if(c>10000){
        break
      }
    }
  }
  return(c)
}




tael1_laplace=function(m=50,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=200){
  false=0
  x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(0.5,times=d)
  y2=rep(-0.5,times=d)
  for(i in 1:m){
    x1[i,]=exp(-y1%*%x[i,])
    x2[i,]=exp(-y2%*%x[i,])
  }
  x3=colMeans(x1)
  x4=colMeans(x2)
  x5=cbind(x3,x4)
  x6=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x6[i,]),exp(-y2%*%x6[i,]))
  }
  x8=colMeans(x7)
  repeat{
    c=d+1
    repeat{
      x9=matrix(0,c,d)
      h0=matrix(0,c,2)
      h1=matrix(0,c,2)
      for(i in 1:c){
        x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
        h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-(d-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          l0=l0$`-2LLR`
          l1=l1$`-2LLR`
          tl0=l0*max(1-l0/nrow(x10),0.5)
          tl1=l1*max(1-l1/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
        }
        if(sum(z)>threshold){
          break
        }
        c=c+1
        if(c>=50){
          break
        }
        x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
      if(sum(z)>threshold){
        false=false+1
        break
      }
      if(c>=50){
        break
      }
    }
    if(c>=50){
      break
    }
    if(false>=1){
      break
    }
  }
  if(false<1){
  nodetect=0
  repeat{
    x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
    h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
    h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    z=numeric()
    for(i in 1:(c-(d-1))){
      x10=h0[i:c,]
      x11=h1[i:c,]
      x10s=colMeans(x10)
      x11s=colMeans(x11)
      x10=rbind(x10,-max(1,log(c)/2)*x10s)
      x11=rbind(x11,-max(1,log(c)/2)*x11s)
      l0=el.test.newton(x=x10,mu=rep(0,2))
      l1=el.test.newton(x=x11,mu=rep(0,2))
      l0=l0$`-2LLR`
      l1=l1$`-2LLR`
      tl0=l0*max(1-l0/nrow(x10),0.5)
      tl1=l1*max(1-l1/nrow(x11),0.5)
      z[i]=exp(0.5*(tl0-tl1))
    }
    c=c+1
    if(sum(z)>threshold){
      break
    }
    if(c>=101){
      break
    }
  }
  if(sum(z)<=threshold){
    repeat{
      x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-(d-1))){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        l0=l0$`-2LLR`
        l1=l1$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      c=c+1
      print(c)
      if(sum(z)>threshold){
        break
      }
      if(c>1000){
        nodetect=nodetect+1
        break
      }
    }
  }
  return(c(c-50,nodetect,false))
  }
  if(false>=1){
    return(c(NA,NA,false))
  }
}


