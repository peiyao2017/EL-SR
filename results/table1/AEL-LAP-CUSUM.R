setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/paper/table 1/ARL1")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)


size=c(1000,200,100,20)
thresholds=c(1.53,1.38,1.14,0.55)


ael1_laplace=function(m=100,mu0=0,mu1=0.5,threshold=20){
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
        for(i in 1:(c-1)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          z[i]=(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
        if(max(z)>threshold){
          break
        }
        c=c+1
        if(c>=50){
          break
        }
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
      if(max(z)>threshold){
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
  if(false==0){
    nodetect=0
    repeat{
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      z=numeric()
      for(i in 1:(c-1)){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        z[i]=(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
      c=c+1
      if(max(z)>threshold){
        break
      }
      if(c>=101){
        break
      }
    }
    if(max(z)<=threshold){
      repeat{
        x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
        h0=h0[-1,]
        h1=h1[-1,]
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-1)){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          z[i]=(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
        print(max(z))
        if(max(z)>threshold){
          break
        }
        c=c+1
        print(c)
        if(max(z)>threshold){
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


ARL1=numeric()
SDRL1=numeric()
ND=numeric()
FAP=numeric()
for(i in 1:length(thresholds)){
  a1=numeric()
  b1=numeric()
  a2=numeric()
  b2=numeric()
  a3=numeric()
  a=foreach(j=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    library(mvtnorm)
    ael1_laplace(threshold=thresholds[i],m=size[i])
  }
  for(j in 1:1000){
    a3[j]=a[3*j]
  }
  repeat{
    a=foreach(j=1:1000, .combine="c")%dopar%{
      library(el.convex)
      library(MASS)
      library(mvtnorm)
      ael1_laplace(threshold=thresholds[i],m=size[i])
    }
    for(j in 1:1000){
      b1[j]=a[3*(j-1)+1]
      b2[j]=a[3*(j-1)+2]
    }
    a1=c(a1,na.omit(b1))
    a2=c(a2,na.omit(b2))
    if(length(a1)>1000&length(a2)>1000){
      break
    }
  }
  a1=a1[1:1000 ]
  a2=a2[1:1000 ]
  ARL1[i]=mean(na.omit(a1))
  SDRL1[i]=sd(na.omit(a1))/sqrt(1000)
  ND[i]=sum(na.omit(a2))
  FAP[i]=sum(na.omit(a3))/1000
}



outname="table1_AEL-LAP-CUSUM.RData"
result=list(data.frame(size,thresholds,ARL1,SDRL1,ND,FAP))
save(result,file = outname)

