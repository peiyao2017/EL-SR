install.packages('doParallel')
library(doParallel)

myCluster <- makeCluster(2, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)



library(el.convex)
ael_mean=function(m=50,mu0=0,mu1=0.5,threshold=1.21,trunction=1000){
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
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,x9[c]-x5)
      h1=rbind(h1,x9[c]-x8)
      if(c>trunction){
        break
      }
    }
  }
  return(c)
}

ael1_mean=function(m=50,mu0=0,mu1=0.5,threshold=2,v=0){
 if(v<=d+1){
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
         if(i <v){
            x9[i]=rnorm(n=1,mean=mu0,sd=1)
         }
         if(i>=v){
           x9[i]=rnorm(n=1,mean=mu1,sd=1)
         }
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
         x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
         h0=c(h0,x9[c]-x5)
         h1=c(h1,x9[c]-x8)
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
          if(c>=v){
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
