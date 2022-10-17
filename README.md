title: "README"
author: "Peiyao"
date: '2022-10-16'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Background

This is an insruction about the code for empirical likelihood based Shiryaev-Roberts change point detection procedure.


## Usage
Before running the functions, please install the packages
```{r}
install.packages("ks")
install.packages("el.convex")
install.packages("MASS")
install.packages("mvtnorm")
library(ks)
library(el.convex)
library(mvtnorm)
library(MASS)


```


For each table in the paper, we use the R functions defined in R files to compute average run lengths. Take table 1  for example, to recover the result in table 1, readers need to change  threshold and sample size $m$ arguements, and leave other arguments unchanged. The same for other tables.\par

```{r}
###Functions for $ARL_0$
#### optimal SR procedure

opt_sr=function(threshold=1,mu0=0,mu1=0.5,covar=1){
  c=0
  r=0
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
    r=(1+r)*dnorm(x,mean=mu1,sd=covar)/dnorm(x,mean=mu0,sd=covar)
    if(r>threshold){
      break
    }
    if(c>1000){
      break
    }
  }
  
  return(c)
}

#### approached SR procedure

app_sr=function(m=100,threshold=100,mu0=0,mu1=0.5,covar=1){
  c=0
  r=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
     r=(1+r)*dnorm(x,mean=estmu1,sd=estcovar1)/dnorm(x,mean=estmu0,sd=estcovar0)
     print(r)
  if(r>threshold){
    break
  }
     if(c>1000){
       break
     }
  }
  return(c)
}

#### normal kernel SR procedure

ks_sr=function(m=100,threshold=1,mu0=0,mu1=0.5,covar=1){
  c=0
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
  r=0
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
      r=(1+r)*f1(x)/f0(x)
      if(r>threshold){
        break
      }
      if(c>1000){
        break
      }
  }
  return(c)
}

#### empirical likelihood SR procedures based on mean

#### orginal EL

el_mean=function(m=50,mu0=0,mu1=0.5,threshold=200){
     x=rnorm(n = m, mean=mu0,sd=1)
     x1=rnorm(n=m,mean=mu1,sd=1)
     x5=mean(x)
     x8=mean(x1)
     c=3
       x9=numeric()
       h0=numeric()
       h1=numeric()
       for(i in 1:c){
           x9[i]=rnorm(n=1,mean=mu0,sd=1)
         }
       for(i in 1:c){
           h0[i]=x9[i]-x5
           h1[i]=x9[i]-x8
         }
       repeat{
           z=numeric()
           for(i in 1:(c-1)){
              x10=h0[i:c]
               x11=h1[i:c]
               l0=el.test.newton(x=x10,mu=rep(0,1))
               l1=el.test.newton(x=x11,mu=rep(0,1))
               if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
                 z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
                 }
             }
            z=na.omit(z)
           if(length(z)==0){
             x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
             c=c+1
             print(c)
             h0=c(h0,x9[c]-x5)
             h1=c(h1,x9[c]-x8)
             }
           if(length(z)>0){
             print(sum(z))
             if(sum(z)<=threshold){
                 x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
                 c=c+1
                 print(c)
                 h0=c(h0,x9[c]-x5)
                 h1=c(h1,x9[c]-x8)
             }
             }
           if(length(z)>0){
             print(sum(z))
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
             z=numeric()
             current0=h0[(length(h0)-99):length(h0)]
             current1=h1[(length(h1)-99):length(h1)]
             for(i in 1:(length(current0)-1)){
                 x10=current0[i:length(current0)]
                 x11=current1[i:length(current0)]
                 l0=el.test.newton(x=x10,mu=rep(0,1))
                 l1=el.test.newton(x=x11,mu=rep(0,1))
                
                 if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
                   z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
                   }
               }
               z=na.omit(z)
               if(length(z)==0){
                   x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
                   c=c+1
                   print(c)
                   h0=c(h0,x9[c]-x5)
                   h1=c(h1,x9[c]-x8)
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
                     x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
                     c=c+1
                     print(c)
                     print(sum(z))
                     h0=c(h0,x9[c]-x5)
                     h1=c(h1,x9[c]-x8)
                 }
               }
               if(c>1000){
                 break
               }
               }
         }
     return(c)
   
     }

#### adjusted EL

ael_mean=function(m=50,mu0=0,mu1=0.5,threshold=1.21){
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rnorm(n = m,mean=mu1,sd=1)
  x5=mean(x)
  x8=mean(x1)
  c=3
  if(c<=100){
    x9=numeric()
    h0=numeric()
    h1=numeric()
    for(i in 1:c){
      x9[i]=rnorm(n=1,mean=mu0,sd=1)
    }
    for(i in 1:c){
      h0[i]=x9[i]-x5
      h1[i]=x9[i]-x8
    }
    z=numeric()
    repeat{
      for(i in 1:c){
        x10=h0[i:c]
        x11=h1[i:c]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0) 
        l1=el.test.newton(x=x11,mu=0) 
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,x9[c]-x5)
      h1=rbind(h1,x9[c]-x8)
      if(c>100){
        break
      }
    }
  }
  if(sum(z)<=threshold){
    repeat{
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      z=numeric()
      for(i in 1:length(current0)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0) 
        l1=el.test.newton(x=x11,mu=0) 
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,x9[c]-x5)
      h1=rbind(h1,x9[c]-x8)
      if(c>1000){
        break
      }
    }
  }
  return(c)
}

#### transformed EL

tel_mean=function(m=50,mu0=0,mu1=0.5,threshold=10){
  x=rnorm(n = m, mean=mu0, sd=1)
  x1=rnorm(n = m, mean=mu1, sd=1)
  x5=mean(x)
  x8=mean(x1)
  c=3
  x9=numeric()
  h0=numeric()
  h1=numeric()
  for(i in 1:c){
    x9[i]=rnorm(n=1,mean=mu0,sd=1)
  }
  for(i in 1:c){
    h0[i]=x9[i]-x5
    h1[i]=x9[i]-x8
  }
  repeat{
    z=numeric()
    for(i in 1:(c-2)){
      x10=h0[i:c]
      x11=h1[i:c]
      l0=el.test.newton(x=x10,mu=0) 
      l1=el.test.newton(x=x11,mu=0) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/length(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/length(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
    }
    if(length(z)>0){
      if(sum(z)<=threshold){
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
      z=numeric()
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      for(i in 1:(length(current0)-2)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        l0=el.test.newton(x=x10,mu=0) 
        l1=el.test.newton(x=x11,mu=0) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/length(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/length(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
      }
      if(length(z)>0){
        if(sum(z)>threshold){
          print(sum(z))
          break
        }
      }
      if(length(z)>0){
        if(sum(z)<=threshold){
          x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
          c=c+1
          print(c)
          print(sum(z))
          h0=c(h0,x9[c]-x5)
          h1=c(h1,x9[c]-x8)
        }
      }
      if(c>1000){
        break
      }
    }
  }
  return(c)
  
}

#### transformed adjusted EL

tael_mean=function(m=50,mu0=0,mu1=0.5,threshold=100){
  x=rnorm(n = m, mean=mu0, sd=1)
  x1=rnorm(n = m, mean=mu1, sd=1)
  x5=mean(x)
  x8=mean(x1)
  c=3
  if(c<=100){
    x9=numeric()
    h0=numeric()
    h1=numeric()
    for(i in 1:c){
      x9[i]=rnorm(n=1,mean=mu0,sd=1)
    }
    for(i in 1:c){
      h0[i]=x9[i]-x5
      h1[i]=x9[i]-x8
    }
    z=numeric()
    repeat{
      for(i in 1:(c-1)){
        x10=h0[i:c]
        x11=h1[i:c]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0)$`-2LLR`
        l1=el.test.newton(x=x11,mu=0)$`-2LLR`
        tl0=l0*max(1-l0/length(x10),0.5)
        tl1=l1*max(1-l1/length(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
      if(c>100){
        break
      }
    }
  }
  if(sum(z)<=threshold){
    repeat{
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      z=numeric()
      for(i in 1:(length(current0)-1)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0)$`-2LLR`
        l1=el.test.newton(x=x11,mu=0)$`-2LLR`
        tl0=l0*max(1-l0/length(x10),0.5)
        tl1=l1*max(1-l1/length(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
      if(c>1000){
        break
      }
    }
  }
  return(c)
}

#### empirical likelihood SR procedures based on Laplace transformation

#### orginal EL

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
      if(c>1000){
        break
      }
    }
  }
  return(c)
}

#### adjusted EL

ael_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
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
   if(c<=100){
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
         z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
       }
       if(sum(z)>threshold){
         break
       }
       x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
       c=c+1
       print(c)
       h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
       h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
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
       for(i in 1:(nrow(current0)-1)){
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
         z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
       }
       print(sum(z))
       if(sum(z)>threshold){
         break
       }
       x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
       c=c+1
       print(c)
       h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
       h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
       h0=h0[-1,]
       h1=h1[-1,]
       if(c>1000){
         break
       }
     }
   }
   return(c)
}

#### transformed EL

tel_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
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
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
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
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
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
      if(c>1000){
        break
      }
    }
  }
  return(c)
}

#### transformed adjusted EL

tael_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
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
  if(c<=100){
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
        if(i==c){
          x10s=x10
          x11s=x11
        }else{
          x10s=colMeans(x10)
          x11s=colMeans(x11)
        }
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
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
      for(i in 1:(nrow(current0)-1)){
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
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
      print(sum(z))
      if(sum(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      if(c>1000){
        break
      }
    }
    
  }
  return(c)
}

### functions to compute $ARL_1$, name rules are the same as above.
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


el1_mean=function(m=100,mu0=0,mu1=0.5,threshold=200){
  false=0
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rnorm(n = m,mean=mu1,sd=1)
  x5=mean(x)
  x8=mean(x1)
  repeat{
    c=3
    repeat{
    x9=numeric()
    h0=numeric()
    h1=numeric()
    for(i in 1:c){
      x9[i]=rnorm(n=1,mean=mu0,sd=1)
    }
    for(i in 1:c){
      h0[i]=x9[i]-x5
      h1[i]=x9[i]-x8
    }
    repeat{
      z=numeric()
      for(i in 1:(c-1)){
        x10=h0[i:c]
        x11=h1[i:c]
        l0=el.test.newton(x=x10,mu=rep(0,1))
        l1=el.test.newton(x=x11,mu=rep(0,1))
         
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
      }
      if(length(z)>0){
        if(sum(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        if(sum(z)<=threshold){
          x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
          c=c+1
          print(c)
          h0=c(h0,x9[c]-x5)
          h1=c(h1,x9[c]-x8)
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
    if(false>=10){
      break
    }
  }
  if(false<10){
  nodetect=0
  repeat{
    z=numeric()
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    c=c+1
    print(c)
    h0=c(h0,x9[c]-x5)
    h1=c(h1,x9[c]-x8)
    z=numeric()
    for(i in 1:(c-1)){
        x10=h0[i:c]
        x11=h1[i:c]
        l0=el.test.newton(x=x10,mu=1)
        l1=el.test.newton(x=x11,mu=1)
        
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
    }
    if(length(z)>0){
      if(sum(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(sum(z)<=threshold){
        x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    c=c+1
    print(c)
    h0=c(h0,x9[c]-x5)
    h1=c(h1,x9[c]-x8)
    current0=h0[(length(h0)-99):length(h0)]
    current1=h1[(length(h1)-99):length(h1)]
    z=numeric()
    for(i in 1:(length(current0)-1)){
      x10=current0[i:length(current0)]
      x11=current1[i:length(current0)]
      l0=el.test.newton(x=x10,mu=0)
      l1=el.test.newton(x=x11,mu=0)
       
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
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
        x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
      }
    }
    if(c>1000){
      nodetect=nodetect+1
      break
    }
  }
  }
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}


ael1_mean=function(m=50,mu0=0,mu1=0.5,threshold=2){
  false=0
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rnorm(n = m, mean=mu1,sd=1)
  x5=mean(x)
  x8=mean(x1)
  repeat{
    c=3
    repeat{
      x9=numeric()
      h0=numeric()
      h1=numeric()
      for(i in 1:c){
        x9[i]=rnorm(n=1,mean=mu0,sd=1)
      }
      for(i in 1:c){
        h0[i]=x9[i]-x5
        h1[i]=x9[i]-x8
      }
      z=numeric()
      repeat{
        for(i in 1:c){
          x10=h0[i:c]
          x11=h1[i:c]
          x10s=mean(x10)
          x11s=mean(x11)
          x10=c(x10,-max(1,log(c)/2)*x10s)
          x11=c(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=0)
          l1=el.test.newton(x=x11,mu=0)
          z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
        }
        if(sum(z)>threshold){
          break
        }
        c=c+1
        if(c>=50){
          break
        }
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
    if(false>=10){
      break
    }
  }
  if(false<10){
  nodetect=0
  repeat{
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    h0=rbind(h0,x9[c]-x5)
    h1=rbind(h1,x9[c]-x8)
    z=numeric()
    for(i in 1:c){
      x10=h0[i:c]
      x11=h1[i:c]
      x10s=mean(x10)
      x11s=mean(x11)
      x10=c(x10,-max(1,log(c)/2)*x10s)
      x11=c(x11,-max(1,log(c)/2)*x11s)
      l0=el.test.newton(x=x10,mu=0)
      l1=el.test.newton(x=x11,mu=0)
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
      x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=rbind(h0,x9[c]-x5)
      h1=rbind(h1,x9[c]-x8)
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      z=numeric()
      for(i in 1:length(current0)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0)
        l1=el.test.newton(x=x11,mu=0)
        z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}


tel1_mean=function(m=100,mu0=0,mu1=0.5,threshold=200){
  false=0
  x=rnorm(n = m, mean=mu0, sd=1)
  x1=rnorm(n = m, mean=mu1, sd=1)
  x5=mean(x)
  x8=mean(x1)
  repeat{
    c=3
    repeat{
      x9=numeric()
      h0=numeric()
      h1=numeric()
      for(i in 1:c){
        x9[i]=rnorm(n=1,mean=mu0,sd=1)
      }
      for(i in 1:c){
        h0[i]=x9[i]-x5
        h1[i]=x9[i]-x8
      }
      repeat{
        z=numeric()
        for(i in 1:(c-2)){
          x10=h0[i:c]
          x11=h1[i:c]
          l0=el.test.newton(x=x10,mu=0) 
          l1=el.test.newton(x=x11,mu=0) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/length(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/length(x11),0.5)
            z[i]=exp(0.5*(tl0-tl1))
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
          c=c+1
          print(c)
          h0=c(h0,x9[c]-x5)
          h1=c(h1,x9[c]-x8)
        }
        if(length(z)>0){
          if(sum(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          if(sum(z)<=threshold){
            x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
            c=c+1
            print(c)
            h0=c(h0,x9[c]-x5)
            h1=c(h1,x9[c]-x8)
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
    if(false>=10){
      break
    }
  }
  nodetect=0
  if(false<10){
  repeat{
    z=numeric()
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    c=c+1
    print(c)
    h0=c(h0,x9[c]-x5)
    h1=c(h1,x9[c]-x8)
    z=numeric()
    for(i in 1:(c-2)){
      x10=h0[i:c]
      x11=h1[i:c]
      l0=el.test.newton(x=x10,mu=0) 
      l1=el.test.newton(x=x11,mu=0) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/length(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/length(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
    }
    if(length(z)>0){
      if(sum(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(sum(z)<=threshold){
        x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      c=c+1
      print(c)
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      z=numeric()
      for(i in 1:(length(current0)-2)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        l0=el.test.newton(x=x10,mu=0) 
        l1=el.test.newton(x=x11,mu=0) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/length(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/length(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
        c=c+1
        print(c)
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
          x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
          c=c+1
          print(c)
          h0=c(h0,x9[c]-x5)
          h1=c(h1,x9[c]-x8)
        }
      }
      if(c>1000){
        nodetect=nodetect+1
        break
      }
    }
  }
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}


tael1_mean=function(m=50,mu0=0,mu1=0.5,threshold=100){
  false=0
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rnorm(n = m, mean=mu1, sd=1)
  x5=mean(x)
  x8=mean(x1)
  repeat{
    c=3
    repeat{
      x9=numeric()
      h0=numeric()
      h1=numeric()
      for(i in 1:c){
        x9[i]=rnorm(n=1,mean=mu0,sd=1)
      }
      for(i in 1:c){
        h0[i]=x9[i]-x5
        h1[i]=x9[i]-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-1)){
          x10=h0[i:c]
          x11=h1[i:c]
          x10s=mean(x10)
          x11s=mean(x11)
          x10=c(x10,-max(1,log(c)/2)*x10s)
          x11=c(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=0)$`-2LLR`
          l1=el.test.newton(x=x11,mu=0)$`-2LLR`
          tl0=l0*max(1-l0/length(x10),0.5)
          tl1=l1*max(1-l1/length(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
        }
        if(sum(z)>threshold){
          break
        }
        c=c+1
        if(c>=50){
          break
        }
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        h0=c(h0,x9[c]-x5)
        h1=c(h1,x9[c]-x8)
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
    if(false>=10){
      break
    }
  }
  if(false<10){
  nodetect=0
  
  repeat{
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    h0=c(h0,x9[c]-x5)
    h1=c(h1,x9[c]-x8)
    z=numeric()
    for(i in 1:(c-1)){
      x10=h0[i:c]
      x11=h1[i:c]
      x10s=mean(x10)
      x11s=mean(x11)
      x10=c(x10,-max(1,log(c)/2)*x10s)
      x11=c(x11,-max(1,log(c)/2)*x11s)
      l0=el.test.newton(x=x10,mu=0)$`-2LLR`
      l1=el.test.newton(x=x11,mu=0)$`-2LLR`
      tl0=l0*max(1-l0/length(x10),0.5)
      tl1=l1*max(1-l1/length(x11),0.5)
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
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=c(h0,x9[c]-x5)
      h1=c(h1,x9[c]-x8)
      current0=h0[(length(h0)-99):length(h0)]
      current1=h1[(length(h1)-99):length(h1)]
      z=numeric()
      for(i in 1:(length(current0)-1)){
        x10=current0[i:length(current0)]
        x11=current1[i:length(current0)]
        x10s=mean(x10)
        x11s=mean(x11)
        x10=c(x10,-max(1,log(c)/2)*x10s)
        x11=c(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=0)$`-2LLR`
        l1=el.test.newton(x=x11,mu=0)$`-2LLR`
        tl0=l0*max(1-l0/length(x10),0.5)
        tl1=l1*max(1-l1/length(x11),0.5)
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
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
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
      if(sum(z)>threshold){
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
  return(c(c-50,nodetect))
  }
  if(false){
    return("false")
  }
}


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
           z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
         }
         if(sum(z)>threshold){
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
     if(false>=10){
       break
     }
   }
   if(false<10){
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
         z[i]=exp(0.5*(l0$`-2LLR`-l1$`-2LLR`))
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
   return(c(c-50,nodetect))
   }
   if(false>=10){
     return("false")
   }
 }
 

tel1_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
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
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
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
    if(sum(z)>threshold){
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
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=exp(0.5*(tl0-tl1))
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
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
        z[i]=exp(0.5*(tl0-tl1))
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
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
  }


tael1_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
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
          l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
          l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
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
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
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
    if(false>=10){
      break
    }
  }
  if(false<10){
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
      l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
      l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
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
      x9=rbind(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
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
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
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
  return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}





```
## Inplementation
The computation is implemented on Ohio super computer's cluster with 80 cores. To use cluster computing, run the code below.


```{r}
install.packages("doParalle")
library(doParallel)
number_of_cores=detectCores()
myCluster <- makeCluster(number_of_cores-1, # number of cores to use
                         type = "PSOCK") # type of cluster
#notice that we need to leave one core for computer system
registerDoParallel(myCluster)
```

Average run lengths are estimated with for loop. The foreach function let you run the loop with multiple cores, and the ".combine="c"" argument combines the result on each core as a vector. Take the original EL as an example. Notice we need to include the required packages in the foreach loop.\par
```{r}

ARL0=foreach(i=1:1000, .combine="c")%dopar%{
  library(el.convex)
  library(MASS)
  library(mvtnorm)
  el_mean(m=20,threshold=100)
}
mean(ARL0)
sd(ARL0)/sqrt(1000)



```
To compute $ARL_1$, we set the change point at the 50th observation. However, when the training sample is too small, the procedures produce to many false alarms before $n=50$. If a false alarm happens, we discard the current sequence and start a new one. If more than 10 false alarms occur, we stop the procedure and the function returns a "false". We repeat the function until we get 1000 simulated detection delays. Take the original EL for an example. The function el1_mean returns two elements, the delay of detection and a indicator telling if the change point is identified before the 1000th observation. If the change point is identified, the indicator is 0. Otherwise it is 1.
```{r}
arl1=numeric()
repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    el1_mean(threshold =250,m=20)
  }
  arl1=c(arl1,a[a!="false"])
  if(length(arl1)>2000){
    break
  }
}

arl1=arl1[1:2000]
arl1=as.numeric(arl1)
delay=numeric()
nodetect=numeric()
for(i in 1:1000){
  delay[i]=arl1[2*i-1]
  nodetect[i]=arl1[2*i]
}
mean(delay)
sd(delay)/sqrt(1000)
sum(nodetect)

```
## Data

The aircraft data mentioned in the paper is available as a spread sheet file on github. The data gives intervals between failures.
