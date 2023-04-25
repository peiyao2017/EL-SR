setwd("/panfs/jay/groups/20/panwei/wan01299/collider_bias/paper/table 4/ARL1")

library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-10, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)


size=c(1000,200,100,20)
thresholds=c(19.7,21,22.9,28)
tel1_var=function(m=100,d=2,threshold=200){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=rgamma(n=d,shape=10,rate=2)
  }
  
  x5=apply(x,MARGIN =2,var)
  x51=colMeans(x)
  x7=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x7[i,]=rgamma(n=d,shape=5,rate=1)
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
        x9[i,]=rgamma(n=d,shape=10,rate=2)
      }
      for(i in 1:c){
        h0[i,]=(x9[i,]-x51)^2-x5
        h1[i,]=(x9[i,]-x81)^2-x8
      }
      repeat{
        z=numeric()
        for(i in 1:(c-d)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,d)) 
          l1=el.test.newton(x=x11,mu=rep(0,d)) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=exp(0.5*(tl0-tl1))
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rgamma(n=d,shape=10,rate=2))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
        }
        if(length(z)>0){
          if(sum(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          if(sum(z)<=threshold){
            x9=rbind(x9,rgamma(n=d,shape=10,rate=2))
            c=c+1
            print(c)
            h0=rbind(h0,(x9[c,]-x51)^2-x5)
            h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
      x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      z=numeric()
      for(i in 1:(c-d)){
        x10=h0[i:c,]
        x11=h1[i:c,]
        l0=el.test.newton(x=x10,mu=rep(0,d)) 
        l1=el.test.newton(x=x11,mu=rep(0,d)) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
      }
      if(length(z)>0){
        if(sum(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        if(sum(z)<=threshold){
          x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
        x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-d)){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          l0=el.test.newton(x=x10,mu=rep(0,d)) 
          l1=el.test.newton(x=x11,mu=rep(0,d)) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=exp(0.5*(tl0-tl1))
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
            x9=rbind(x9,rgamma(n=d,shape=5,rate=1))
            c=c+1
            print(c)
            h0=rbind(h0,(x9[c,]-x51)^2-x5)
            h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
    tel1_var(threshold=thresholds[i],m=size[i])
  }
  for(j in 1:1000){
    a3[j]=a[3*j]
  }
  repeat{
    a=foreach(j=1:1000, .combine="c")%dopar%{
      library(el.convex)
      library(MASS)
      library(mvtnorm)
      tel1_var(threshold=thresholds[i],m=size[i])
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









outname="table4_TEL-MEAN and VAR.RData"
result=list(data.frame(size,thresholds,ARL1,SDRL1,ND,FAP))
save(result,file = outname)

