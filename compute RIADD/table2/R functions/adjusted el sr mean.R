install.packages('doParallel')
library(doParallel)

myCluster <- makeCluster(2, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)


library(el.convex)


ael_mean=function(m=100,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=100,trunction=1000){
  x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
  x1=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
  x5=colMeans(x)
  x8=colMeans(x1)
  c=d+1
  if(c<=100){
    x9=matrix(0,c,d)
    h0=matrix(0,c,d)
    h1=matrix(0,c,d)
    for(i in 1:c){
      x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
    }
    for(i in 1:c){
      h0[i,]=x9[i,]-x5
      h1[i,]=x9[i,]-x8
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
        l0=el.test.newton(x=x10,mu=rep(0,d)) 
        l1=el.test.newton(x=x11,mu=rep(0,d)) 
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,x9[c,]-x5)
      h1=rbind(h1,x9[c,]-x8)
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
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,d)) 
        l1=el.test.newton(x=x11,mu=rep(0,d)) 
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,x9[c,]-x5)
      h1=rbind(h1,x9[c,]-x8)
      if(c>trunction){
        break
      }
    }
  }
  return(c)
}

ael1_mean=function(m=100,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=100,v=0){
  if(v<=d+1){
    x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
    x1=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
    x5=colMeans(x)
    x8=colMeans(x1)
    repeat{
      c=d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,d)
        h1=matrix(0,c,d)
        for(i in 1:c){
          if(i<c){
            x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
          }
          if(i>=c){
            x9[i,]=mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d))
          }
        }
        for(i in 1:c){
          h0[i,]=x9[i,]-x5
          h1[i,]=x9[i,]-x8
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
            l0=el.test.newton(x=x10,mu=rep(0,d))
            l1=el.test.newton(x=x11,mu=rep(0,d))
            z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
          }
          if(sum(z)>threshold){
            break
          }
          c=c+1
          if(c>=50){
            break
          }
          x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
          h0=rbind(h0,x9[c,]-x5)
          h1=rbind(h1,x9[c,]-x8)
        }
        if(sum(z)>threshold){
           
          break
        }
        if(c>=50){
          break
        }
      }
      if(c>=50){
        break
      }
      if(sum(z)>threshold){
        break
      }
    }
    if(sum(z)>threshold){
      delay=c-v
      false=0
      nodetect=0
      return(c(delay,false,nodetect))
    }
    
      nodetect=0
      repeat{
        x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
        h0=rbind(h0,x9[c,]-x5)
        h1=rbind(h1,x9[c,]-x8)
        z=numeric()
        for(i in 1:(c-(d-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,d))
          l1=el.test.newton(x=x11,mu=rep(0,d))
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
          h0=rbind(h0,x9[c,]-x5)
          h1=rbind(h1,x9[c,]-x8)
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
            l0=el.test.newton(x=x10,mu=rep(0,d))
            l1=el.test.newton(x=x11,mu=rep(0,d))
            z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
  if(v>d+1){
    false=0
    x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
    x1=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
    x5=colMeans(x)
    x8=colMeans(x1)
    repeat{
      c=d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,d)
        h1=matrix(0,c,d)
        for(i in 1:c){
          x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
        }
        for(i in 1:c){
          h0[i,]=x9[i,]-x5
          h1[i,]=x9[i,]-x8
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
            l0=el.test.newton(x=x10,mu=rep(0,d))
            l1=el.test.newton(x=x11,mu=rep(0,d))
            z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
          }
          if(sum(z)>threshold){
            break
          }
          c=c+1
          if(c>=v){
            break
          }
          x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
          h0=rbind(h0,x9[c,]-x5)
          h1=rbind(h1,x9[c,]-x8)
        }
        if(sum(z)>threshold){
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
      if(false>=1){
        break
      }
    }
    if(false<1){
      nodetect=0
      repeat{
        x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
        h0=rbind(h0,x9[c,]-x5)
        h1=rbind(h1,x9[c,]-x8)
        z=numeric()
        for(i in 1:(c-(d-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,d))
          l1=el.test.newton(x=x11,mu=rep(0,d))
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
          h0=rbind(h0,x9[c,]-x5)
          h1=rbind(h1,x9[c,]-x8)
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
            l0=el.test.newton(x=x10,mu=rep(0,d))
            l1=el.test.newton(x=x11,mu=rep(0,d))
            z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
      return(c(c-v,false,nodetect))
    }
    if(false>=1){
      return(c(NA,false,NA))
    }
  }  
}
