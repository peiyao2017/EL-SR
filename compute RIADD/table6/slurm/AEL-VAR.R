setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/paper/computeRAIDD/table6/")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)


vs=c(0:20)
size=c(1000,200,100,20,20)
thresholds=c(43.5,57.5,60,58,2350)
N=c(100,100,100,100,1000)
trunctions=c(1000,1000,1000,1000,3000)
ael_mean=function(m=100,d=2,lambda0=4,lambda1=4.5,threshold=200,trunction=1000){
  x=matrix(0,nrow=m,ncol=d)
  y=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=rpois(n=d,lambda=lambda0)
    y[i,]=rpois(n=d,lambda=lambda1)
  }
  x5=colMeans(x)
  x8=colMeans(y)
  x51=numeric()
  x81=numeric()
  for(i in 1:d){
    x51[i]=var(x[,i])
    x81[i]=var(y[,i])
  }
  c=2*d+1
  if(c<=100){
    x9=matrix(0,c,d)
    h0=matrix(0,c,2*d)
    h1=matrix(0,c,2*d)
    for(i in 1:c){
      x9[i,]=rpois(n=d,lambda=lambda0)
    }
    for(i in 1:c){
      h0[i,]=c(x9[i,]-x5,(x9[i,]-x5)^2-x51)
      h1[i,]=c(x9[i,]-x8,(x9[i,]-x8)^2-x81)
    }
    z=numeric()
    repeat{
      for(i in 1:(c-(2*d-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2*d)) 
        l1=el.test.newton(x=x11,mu=rep(0,2*d)) 
        z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,rpois(n=d,lambda=lambda0))
      c=c+1
      print(c)
      h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
      h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
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
      for(i in 1:(nrow(current0)-(2*d-1))){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2*d)) 
        l1=el.test.newton(x=x11,mu=rep(0,2*d)) 
        z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,rpois(n=d,lambda=lambda0))
      c=c+1
      print(c)
      h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
      h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
      if(c>trunction){
        break
      }
    }
  }
  return(c)
}

ael1_mean=function(m=100,d=2,lambda0=4,lambda1=4.5,threshold=200,v=0){
  if(v<=2*d+1){
    x=matrix(0,ncol = d,nrow=m)
    y=matrix(0,ncol = d,nrow=m)
    for(i in 1:m){
      x[i,]=rpois(n=d,lambda=lambda0)
      y[i,]=rpois(n=d,lambda=lambda1)
    }
    x5=colMeans(x)
    x8=colMeans(y)
    x51=numeric()
    x81=numeric()
    for(i in 1:d){
      x51[i]=var(x[,i])
      x81[i]=var(y[,i])
    }
    repeat{
      c=2*d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,2*d)
        h1=matrix(0,c,2*d)
        for(i in 1:c){
          if(c<v){
            x9[i,]=rpois(n=d,lambda=lambda0)
          }
          if(c>=v){
            x9[i,]=rpois(n=d,lambda=lambda1)
          }
        }
        for(i in 1:c){
          h0[i,]=c(x9[i,]-x5,(x9[i,]-x5)^2-x51)
          h1[i,]=c(x9[i,]-x8,(x9[i,]-x8)^2-x81)
        }
        z=numeric()
        repeat{
          for(i in 1:(c-(2*d-1))){
            x10=h0[i:c,]
            x11=h1[i:c,]
            x10s=colMeans(x10)
            x11s=colMeans(x11)
            x10=rbind(x10,-max(1,log(c)/2)*x10s)
            x11=rbind(x11,-max(1,log(c)/2)*x11s)
            l0=el.test.newton(x=x10,mu=rep(0,2*d))
            l1=el.test.newton(x=x11,mu=rep(0,2*d))
            z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          }
          if(sum(z)>threshold){
            break
          }
          c=c+1
          if(c>=50){
            break
          }
          x9=rbind(x9,rpois(n=d,lambda=lambda1))
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
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
    repeat{
      x9=rbind(x9,rpois(n=d,lambda=lambda1))
      h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
      h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
      z=numeric()
      for(i in 1:(c-(2*d-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2*d))
        l1=el.test.newton(x=x11,mu=rep(0,2*d))
        z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
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
        x9=rbind(x9,rpois(n=d,lambda=lambda1))
        h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
        h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-(2*d-1))){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2*d))
          l1=el.test.newton(x=x11,mu=rep(0,2*d))
          z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
        }
        print(sum(z))
        if(sum(z)>threshold){
          break
        }
        c=c+1
        print(c)
        if(sum(z)>threshold){
          break
        }
        if(c>1000){
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
  if(v>2*d+1){
    false=0
    x=matrix(0,ncol = d,nrow=m)
    y=matrix(0,ncol = d,nrow=m)
    for(i in 1:m){
      x[i,]=rpois(n=d,lambda=lambda0)
      y[i,]=rpois(n=d,lambda=lambda1)
    }
    x5=colMeans(x)
    x8=colMeans(y)
    x51=numeric()
    x81=numeric()
    for(i in 1:d){
      x51[i]=var(x[,i])
      x81[i]=var(y[,i])
    }
    repeat{
      c=2*d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,2*d)
        h1=matrix(0,c,2*d)
        for(i in 1:c){
          x9[i,]=rpois(n=d,lambda=lambda0)
        }
        for(i in 1:c){
          h0[i,]=c(x9[i,]-x5,(x9[i,]-x5)^2-x51)
          h1[i,]=c(x9[i,]-x8,(x9[i,]-x8)^2-x81)
        }
        z=numeric()
        repeat{
          for(i in 1:(c-(2*d-1))){
            x10=h0[i:c,]
            x11=h1[i:c,]
            x10s=colMeans(x10)
            x11s=colMeans(x11)
            x10=rbind(x10,-max(1,log(c)/2)*x10s)
            x11=rbind(x11,-max(1,log(c)/2)*x11s)
            l0=el.test.newton(x=x10,mu=rep(0,2*d))
            l1=el.test.newton(x=x11,mu=rep(0,2*d))
            z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          }
          if(sum(z)>threshold){
            false=false+1
            break
          }
          c=c+1
          if(c>=v){
            break
          }
          x9=rbind(x9,rpois(n=d,lambda=lambda0))
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        }
        if(false>=1){
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
        x9=rbind(x9,rpois(n=d,lambda=lambda1))
        h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
        h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        z=numeric()
        for(i in 1:(c-(2*d-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2*d))
          l1=el.test.newton(x=x11,mu=rep(0,2*d))
          z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
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
          x9=rbind(x9,rpois(n=d,lambda=lambda1))
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
          current0=h0[(nrow(h0)-99):nrow(h0),]
          current1=h1[(nrow(h1)-99):nrow(h1),]
          z=numeric()
          for(i in 1:(nrow(current0)-(2*d-1))){
            x10=current0[i:nrow(current0),]
            x11=current1[i:nrow(current0),]
            x10s=colMeans(x10)
            x11s=colMeans(x11)
            x10=rbind(x10,-max(1,log(c)/2)*x10s)
            x11=rbind(x11,-max(1,log(c)/2)*x11s)
            l0=el.test.newton(x=x10,mu=rep(0,2*d))
            l1=el.test.newton(x=x11,mu=rep(0,2*d))
            z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          }
          print(sum(z))
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
      return(c(c-v,false,nodetect ))
    }
    if(false>=1){
      return(c(NA,false,NA ))
    }
  }
}

sddelay=numeric()
delay=numeric()
for(i in 1:length(size)){
  c1=foreach(l=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    library(mvtnorm)
    ael_mean(threshold=thresholds[i],m=size[i],trunction=trunctions[i] )
  }
  p=numeric()
  sddelay1=numeric()
  b1=numeric()
  for(j in 1:length(vs)){
    a1=numeric()
    a11=numeric()
    
    repeat{
      a=foreach(l=1:1000, .combine="c")%dopar%{
        library(el.convex)
        library(MASS)
        library(mvtnorm)
        ael1_mean(threshold=thresholds[i],m=size[i],v=vs[j])
      }
      for(k in 1:1000){
        a11[k]=a[3*(k-1)+1]
      }
      a11=na.omit(a11)
      a1=c(a1,a11)
      if(length(a1)>1000){
        break
      }
    }
    a1=a1[1:1000]
    b1[j]=mean(a1)
  
    
    p[j]=sum(c1>vs[j])/length(c1)
    p1=as.numeric(c(c1>vs[j]))
    sddelay1[j]=sqrt(var(a1)*mean(p1^2)+mean(a1)^2*var(p1))/sqrt(1000*N[i])
  }
  delay[i]=sum(b1*p)/N[i]
  sddelay[i]=sum(sddelay1)
}

result=list(data.frame(size,thresholds,delay,sddelay))
save(result,file="AEL-VAR.RData")
