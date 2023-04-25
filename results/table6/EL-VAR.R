setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/paper/table 6/ARL1")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)


size=c(1000,200,100,20,20)
thresholds=c(4000,7000,15000,270000,1.2e+23)

el1_mean=function(m=100,d=2,lambda0=4,lambda1=4.5,threshold=200){
  false=0
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
      repeat{
        z=numeric()
        for(i in 1:(c-2*d)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,2*d))
          l1=el.test.newton(x=x11,mu=rep(0,2*d))
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rpois(n=d,lambda=lambda0))
          c=c+1
          print(c)
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        }
        if(length(z)>0){
          if(sum(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          if(sum(z)<=threshold){
            x9=rbind(x9,rpois(n=d,lambda=lambda0))
            c=c+1
            print(c)
            h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
            h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
          }
        }
        if(c>=49){
          break
        }
      }
      if(length(z)>0){
        if(sum(z)>threshold){
          false=false+1
          break
        }
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
  if(false<1){
    nodetect=0
    repeat{
      z=numeric()
      x9=rbind(x9,rpois(n=d,lambda=lambda1))
      c=c+1
      print(c)
      h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
      h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
      z=numeric()
      for(i in 1:(c-2*d)){
        x10=h0[i:c,]
        x11=h1[i:c,]
        l0=el.test.newton(x=x10,mu=rep(0,2*d))
        l1=el.test.newton(x=x11,mu=rep(0,2*d))
        
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,rpois(n=d,lambda=lambda1))
        c=c+1
        print(c)
        h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
        h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
      }
      if(length(z)>0){
        if(sum(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        if(sum(z)<=threshold){
          x9=rbind(x9,rpois(n=d,lambda=lambda1))
          c=c+1
          print(c)
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
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
        x9=rbind(x9,rpois(n=d,lambda=lambda1))
        c=c+1
        print(c)
        h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
        h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-2*d)){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          l0=el.test.newton(x=x10,mu=rep(0,2*d))
          l1=el.test.newton(x=x11,mu=rep(0,2*d))
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            z[i]=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rpois(n=d,lambda=lambda1))
          c=c+1
          print(c)
          h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
          h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
        }
        if(length(z)>0){
          if(sum(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          print(sum(z))
          if(sum(z)<=threshold){
            x9=rbind(x9,rpois(n=d,lambda=lambda1))
            c=c+1
            print(c)
            h0=rbind(h0,c(x9[c,]-x5,(x9[c,]-x5)^2-x51))
            h1=rbind(h1,c(x9[c,]-x8,(x9[c,]-x8)^2-x81))
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
    el1_mean(threshold=thresholds[i],m=size[i])
  }
  for(j in 1:1000){
    a3[j]=a[3*j]
  }
  repeat{
    a=foreach(j=1:1000, .combine="c")%dopar%{
      library(el.convex)
      library(MASS)
      library(mvtnorm)
      el1_mean(threshold=thresholds[i],m=size[i])
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









outname="table6_EL-VAR.RData"
result=list(data.frame(size,thresholds,ARL1,SDRL1,ND,FAP))
save(result,file = outname)

