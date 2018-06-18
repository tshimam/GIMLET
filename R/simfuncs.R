sim1 <- function(n,type="linear",model=1,snr1=1/3,snr2=3){
  x <- runif(n)
  u <- runif(n)
  if(type=="linear"){
    signal <- x
  }
  if(type=="parabolic"){
    signal <- 4 * (x - 0.5)^2
  }
  if(type=="cubic"){
    signal <- 128 * (x-1/3)^3 - 48 * (x-1/3)^2 - 12 * (x-1/3) + 2
  }
  if(type=="exponential"){
    signal <- 2^(10 * x) - 1
  }
  if(type=="sinusoidal"){
    signal <- sin(2*pi*x*(1+x))
  }
  if(type=="categorical"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 1/3){
        signal[i] <- 2/3
      } else if (x[i] < 2/3){
        signal[i] <- 1/3
      } else {
        signal[i] <- 1
      }
    }
  }
  if(type=="circle"){
    signal <- (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))
  }
  if(type=="atan"){
    signal <- atan(10*x)
  }
  if(type=="sigmoid"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 49/100){
        signal[i] <- 0
      } else if (x[i] < 51/100){
        signal[i] <- 50*(x[i] - 1/2) + 1/2
      } else {
        signal[i] <- 1
      }
    }
  }
  if(type=="spike"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 1/20){
        signal[i] <- 20 * x[i]
      } else if (x[i] < 1/10){
        signal[i] <- -18 * x[i] + 19 / 100
      } else {
        signal[i] <- - x[i] / 9 + 1 / 9
      }
    }
  }
  if(type=="random"){
    signal <- rnorm(n)
  }
  noise1 <- rnorm(n, sd = snr1 * sd(signal))
  if(model==1){
    noise2 <- (1 - u) * rnorm(n) * snr2 * sd(signal)
  }
  if(model==2){
    noise2 <- runif(n) * rnorm(n) * snr2 * sd(signal)
  }
  noise3 <- rnorm(n, sd = snr1 * sd(u))
  y <- signal + noise1 + noise2
  z <- u + noise3
  return(list(x=x,y=y,z=z,signal=signal,u=u,noise1=noise1,noise2=noise2,noise3=noise3))
}

sim2 <- function(n,type="linear",model=1,snr1=1/3,snr2=3){
  x <- runif(n)
  u1 <- runif(n)
  u2 <- runif(n)
  if(type=="linear"){
    signal <- x
  }
  if(type=="parabolic"){
    signal <- 4 * (x - 0.5)^2
  }
  if(type=="cubic"){
    signal <- 128 * (x-1/3)^3 - 48 * (x-1/3)^2 - 12 * (x-1/3) + 2
  }
  if(type=="exponential"){
    signal <- 2^(10 * x) - 1
  }
  if(type=="sinusoidal"){
    signal <- sin(2*pi*x*(1+x))
  }
  if(type=="categorical"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 1/5){
        signal[i] <- 0.4
      } else if (x[i] < 2/5){
        signal[i] <- 0.8
      } else if (x[i] < 3/5){
        signal[i] <- 0.2
      } else if (x[i] < 4/5){
        signal[i] <- 1.0
      } else {
        signal[i] <- 0.6
      }
    }
  }
  if(type=="circle"){
    signal <- (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))
  }
  if(type=="atan"){
    signal <- atan(10*x)
  }
  if(type=="sigmoid"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 49/100){
        signal[i] <- 0
      } else if (x[i] < 51/100){
        signal[i] <- 50*(x[i] - 1/2) + 1/2
      } else {
        signal[i] <- 1
      }
    }
  }
  if(type=="spike"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 1/20){
        signal[i] <- 20 * x[i]
      } else if (x[i] < 1/10){
        signal[i] <- -18 * x[i] + 19 / 100
      } else {
        signal[i] <- - x[i] / 9 + 1 / 9
      }
    }
  }
  if(type=="random"){
    signal <- rnorm(n)
  }
  noise1 <- rnorm(n, sd = snr1 * sd(signal))
  if(model==1){
    noise2 <- (1 - u1) * (1 - u2) * rnorm(n) * snr2 * sd(signal)
  }
  if(model==2){
    noise2 <- rnorm(n) * snr2 * sd(signal)
  }
  noise3 <- rnorm(n, sd = snr1 * sd(u1))
  noise4 <- rnorm(n, sd = snr1 * sd(u2))
  y <- signal + noise1 + noise2
  z1 <- u1 + noise3
  z2 <- u2 + noise4
  z <- cbind(z1,z2)
  return(list(x=x,y=y,z=z,signal=signal,u1=u1,u2=u2,noise1=noise1,noise2=noise2,noise3=noise3,noise4=noise4))
}

sim3 <- function(n,type="linear",model=1,snr1=1/3,snr2=3){
  x1 <- runif(n)
  x2 <- runif(n)
  x <- x1 * x2
  u <- runif(n)
  if(type=="linear"){
    signal <- x
  }
  if(type=="parabolic"){
    signal <- 4 * (x - 0.5)^2
  }
  if(type=="cubic"){
    signal <- 128 * (x-1/3)^3 - 48 * (x-1/3)^2 - 12 * (x-1/3) + 2
  }
  if(type=="exponential"){
    signal <- 2^(10 * x) - 1
  }
  if(type=="sinusoidal"){
    signal <- sin(2*pi*x*(1+x))
  }
  if(type=="categorical"){
    signal <- rep(0,n)
    for(i in 1:n){
      for(j in 1:n){
        if(x1[i] < 1/3 & x2[i] < 1/3){
          signal[i] <- 2/3
        } else if(x1[i] < 1/3 & x2[i] < 2/3){
          signal[i] <- 0
        } else if(x1[i] < 1/3 & x2[i] < 1){
          signal[i] <- 0
        } else if(x1[i] < 2/3 & x2[i] < 1/3){
          signal[i] <- 0
        } else if(x1[i] < 2/3 & x2[i] < 2/3){
          signal[i] <- 1/2
        } else if(x1[i] < 2/3 & x2[i] < 1){
          signal[i] <- 0
        } else if(x1[i] < 1 & x2[i] < 1/3){
          signal[i] <- 0
        } else if(x1[i] < 1 & x2[i] < 2/3){
          signal[i] <- 0
        } else {
          signal[i] <- 1
        }
      }
    }
  }
  if(type=="circle"){
    signal <- (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))
  }
  if(type=="atan"){
    signal <- atan(10*x)
  }
  if(type=="sigmoid"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 49/100){
        signal[i] <- 0
      } else if (x[i] < 51/100){
        signal[i] <- 50*(x[i] - 1/2) + 1/2
      } else {
        signal[i] <- 1
      }
    }
  }
  if(type=="spike"){
    signal <- rep(0,n)
    for(i in 1:n){
      if(x[i] < 1/20){
        signal[i] <- 20 * x[i]
      } else if (x[i] < 1/10){
        signal[i] <- -18 * x[i] + 19 / 10
      } else {
        signal[i] <- - x[i] / 9 + 1 / 9
      }
    }
  }
  if(type=="random"){
    signal <- rnorm(n)
  }
  noise1 <- rnorm(n, sd = snr1 * sd(signal))
  if(model==1){
    noise2 <- (1 - u) * rnorm(n) * snr2 * sd(signal)
  }
  if(model==2){
    noise2 <- runif(n) * rnorm(n) * snr2 * sd(signal)
  }
  noise3 <- rnorm(n, sd = snr1 * sd(u))
  y <- signal + noise1 + noise2
  z <- u + noise3
  return(list(x=cbind(x1,x2),y=y,z=z,signal=signal,u=u,noise1=noise1,noise2=noise2,noise3=noise3))
}

display.sim3 <- function(n=20,border="gray",col=ramp.col(c("white", "black"))){
  M <- mesh(seq(from=0,to=1,length=n),seq(from=0,to=1,length=n))
  x1 <- M$x
  x2 <- M$y
  #linear
  x <- x1 * x2
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  signal <- x
  persp3D(z = signal, col = col, border = border,main="linear")
  signal <- 4 * (x - 0.5)^2
  persp3D(z = signal, col = col, border = border,main="parabolic")
  signal <- 128 * (x-1/3)^3 - 48 * (x-1/3)^2 - 12 * (x-1/3) + 2
  persp3D(z = signal, col = col, border = border,main="cubic")
  signal <- 2^(10 * x) - 1
  persp3D(z = signal, col = col, border = border,main="exponential")
  signal <- sin(2*pi*x*(1+x))
  persp3D(z = signal, col = col, border = border,main="sinusoidal")
  #signal <- (2*rbinom(nrow(x),1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))
  #persp3D(z = signal, col = col, border = border,main="circle")
  signal <- atan(10*x)
  persp3D(z = signal, col = col, border = border,main="arc tangent")
  signal <- matrix(0,nrow(x),ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(x1[i,j] < 1/3 & x2[i,j] < 1/3){
        signal[i,j] <- 2/3
      } else if(x1[i,j] < 1/3 & x2[i,j] < 2/3){
        signal[i,j] <- 0
      } else if(x1[i,j] < 1/3 & x2[i,j] < 1){
        signal[i,j] <- 0
      } else if(x1[i,j] < 2/3 & x2[i,j] < 1/3){
        signal[i,j] <- 0
      } else if(x1[i,j] < 2/3 & x2[i,j] < 2/3){
        signal[i,j] <- 1/3
      } else if(x1[i,j] < 2/3 & x2[i,j] < 1){
        signal[i,j] <- 0
      } else if(x1[i,j] < 1 & x2[i,j] < 1/3){
        signal[i,j] <- 0
      } else if(x1[i,j] < 1 & x2[i,j] < 2/3){
        signal[i,j] <- 0
      } else {
        signal[i,j] <- 1
      }
    }
  }
  persp3D(z = signal, col = col, border = border,main="categorical")
  signal <- matrix(0,nrow(x),ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(x[i,j] < 49/100){
        signal[i,j] <- 0
      } else if (x[i,j] < 51/100){
        signal[i,j] <- 50*(x[i,j] - 1/2) + 1/2
      } else {
        signal[i,j] <- 1
      }
    }
  }
  persp3D(z = signal, col = col, border = border,main="sigmoid")
  signal <- matrix(0,nrow(x),ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(x[i,j] < 1/20){
        signal[i,j] <- 20 * x[i,j]
      } else if (x[i,j] < 1/10){
        signal[i,j] <- -18 * x[i,j] + 19 / 10
      } else {
        signal[i,j] <- - x[i,j] / 9 + 1 / 9
      }
    }
  }
  persp3D(z = signal, col = col, border = border,main="spike")
}

display.sim1 <- function(n=10000){
  x <- seq(from=0,to=1,length=n)
  #linear
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  signal <- x
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="linear")
  signal <- 4 * (x - 0.5)^2
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="parabolic")
  signal <- 128 * (x-1/3)^3 - 48 * (x-1/3)^2 - 12 * (x-1/3) + 2
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="cubic")
  signal <- 2^(10 * x) - 1
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="exponential")
  signal <- sin(2*pi*x*(1+x))
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="sinusoidal")
  #signal <- (2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))
  #plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="circle")
  signal <- atan(10*x)
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="arc tangent")
  signal <- rep(0,n)
  for(i in 1:n){
    if(x[i] < 1/3){
      signal[i] <- 2/3
    } else if (x[i] < 2/3){
      signal[i] <- 1/3
    } else {
      signal[i] <- 1
    }
  }
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="categorical")
  signal <- rep(0,n)
  for(i in 1:n){
    if(x[i] < 49/100){
      signal[i] <- 0
    } else if (x[i] < 51/100){
      signal[i] <- 50*(x[i] - 1/2) + 1/2
    } else {
      signal[i] <- 1
    }
  }
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="sigmoid")
  signal <- rep(0,n)
  for(i in 1:n){
    if(x[i] < 1/20){
      signal[i] <- 20 * x[i]
    } else if (x[i] < 1/10){
      signal[i] <- -18 * x[i] + 19 / 10
    } else {
      signal[i] <- - x[i] / 9 + 1 / 9
    }
  }
  plot(x=x,y=signal,type="l",pch=19,cex=0.1,main="spike")
}
