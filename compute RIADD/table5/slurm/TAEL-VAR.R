

tael1_var=function(m=100,d=4,min0=-3,max0=3,min1=-3.5,max1=3.5,threshold=200,v=0){
  if(v<=d+1){
    x=matrix(0,nrow=m,ncol=d)
    for(i in 1:m){
      x[i,]=runif(n=d,min=min0,max=max0)
    }
    x5=apply(x,MARGIN =2,var)
    x51=colMeans(x)
    x7=matrix(0,nrow=m,ncol=d)
    for(i in 1:m){
      x7[i,]=runif(n=d,min=min1,max=max1)
    }
    x81=colMeans(x7)
    x8=apply(x7,MARGIN =2,var)
    repeat{
      c=d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,d)
        h1=matrix(0,c,d)
        for(i in 1:c){
          if(i<v){
            x9[i,]=runif(n=d,min=min0,max=max0)
          }
          if(i>=v){
            x9[i,]=runif(n=d,min=min1,max=max1)
          }
        }
        for(i in 1:c){
          h0[i,]=(x9[i,]-x51)^2-x5
          h1[i,]=(x9[i,]-x81)^2-x8
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
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
          l0=l0$`-2LLR`
          l1=l1$`-2LLR`
          tl0=l0*max(1-l0/nrow(x10),0.5)
          tl1=l1*max(1-l1/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
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
  if(v>d+1){
    false=0
    x=matrix(0,nrow=m,ncol=d)
    for(i in 1:m){
      x[i,]=runif(n=d,min=min0,max=max0)
    }
    x5=apply(x,MARGIN =2,var)
    x51=colMeans(x)
    x7=matrix(0,nrow=m,ncol=d)
    for(i in 1:m){
      x7[i,]=runif(n=d,min=min1,max=max1)
    }
    x81=colMeans(x7)
    x8=apply(x7,MARGIN =2,var)
    repeat{
      c=d+1
      repeat{
        x9=matrix(0,c,d)
        h0=matrix(0,c,d)
        h1=matrix(0,c,d)
        for(i in 1:c){
          x9[i,]=runif(n=d,min=min0,max=max0)
        }
        for(i in 1:c){
          h0[i,]=(x9[i,]-x51)^2-x5
          h1[i,]=(x9[i,]-x81)^2-x8
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
            l0=l0$`-2LLR`
            l1=l1$`-2LLR`
            tl0=l0*max(1-l0/nrow(x10),0.5)
            tl1=l1*max(1-l1/nrow(x11),0.5)
            z[i]=exp(0.5*(tl0-tl1))
          }
          if(sum(z)>threshold){
            false=false+1
            break
          }
          c=c+1
          if(c>=v){
            break
          }
          x9=rbind(x9,runif(n=d,min=min0,max=max0))
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
            l0=l0$`-2LLR`
            l1=l1$`-2LLR`
            tl0=l0*max(1-l0/nrow(x10),0.5)
            tl1=l1*max(1-l1/nrow(x11),0.5)
            z[i]=exp(0.5*(tl0-tl1))
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
      return(c(c-v,false,nodetect))
    }
    if(false>=1){
      return(c(NA,false,NA))
    }
  }
}



tael_var=function(m=100,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5,trunction=1000){
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  
  x5=apply(x,MARGIN =2,var)
  x51=colMeans(x)
  x7=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x7[i,]=runif(n=d,min=min1,max=max1)
  }
  x81=colMeans(x7)
  x8=apply(x7,MARGIN =2,var)
  c=d+1
  if(c<=100){
    x9=matrix(0,c,d)
    h0=matrix(0,c,d)
    h1=matrix(0,c,d)
    for(i in 1:c){
      x9[i,]=runif(n=d,min=min0,max=max0)
    }
    for(i in 1:c){
      h0[i,]=(x9[i,]-x51)^2-x5
      h1[i,]=(x9[i,]-x81)^2-x8
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
        l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,runif(n=d,min=min0,max=max0))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
        l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=rbind(x9,runif(n=d,min=min0,max=max0))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      if(c>=trunction){
        break
      }
    }
  }
  return(c)
}




setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/paper/computeRAIDD/table6/")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)


vs=c(0:20)
size=c(1000,200,100,20,20)
thresholds=c(26.5,28 ,29.5,38.9,92.4)
N=c(100,100,100,100,1000)
trunctions=c(1000,1000,1000,1000,3000)

delay=numeric()
sddelay=numeric()
for(i in 1:length(size)){
  c1=foreach(l=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    library(mvtnorm)
    tael_var(threshold=thresholds[i],m=size[i],trunction=trunctions[i] )
  }
  b1=numeric()
  p=numeric()
  sddelay1=numeric()
  for(j in 1:length(vs)){
    a1=numeric()
    a11=numeric()
    
    repeat{
      a=foreach(l=1:1000, .combine="c")%dopar%{
        library(el.convex)
        library(MASS)
        library(mvtnorm)
        tael1_var(threshold=thresholds[i],m=size[i],v=vs[j])
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
save(result,file="TAEL-VAR.RData")
