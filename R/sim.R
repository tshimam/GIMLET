sim <- function(n,sim.typ=1,fun.typ=1,model=1,r1=1,r2=3,
signal.thres=0.25,noise.thres=0.75){
  refx <- seq(from=0,to=1,length=1000)
  if(sim.typ==1){
    x <- runif(n)
    z <- runif(n)
  }
  if(sim.typ==2){
    x <- runif(n)
    z1 <- runif(n)
    z2 <- runif(n)
    z <- z1*z2
  }
  if(sim.typ==3){
    x1 <- runif(n)
    x2 <- runif(n)
    z <- runif(n)
    x <- x1*x2
  }
  noise.level <- rep(0,n)
  if(model==1){
    u <- z
  } else {
    u <- runif(z)
  }
  for(i in 1:n){
    if(u[i] > noise.thres){
      noise.level[i] <- r1
    } else {
      noise.level[i] <- r2
    }
  }
  if(fun.typ==1){
    signalfun <- function(x){x}
  }  #parabolic+noise
  if(fun.typ==2){
    signalfun <- function(x){4*(x-.5)^2}
  }
  if(fun.typ==3){
    signalfun <- function(x){128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)}
  }
  if(fun.typ==4){
    signalfun <- function(x){sin(4*pi*x)}
  }
  if(fun.typ==5){
    signalfun <- function(x){sin(16*pi*x)}
  }
  if(fun.typ==6){
    signalfun <- function(x){x^(1/4)}
  }
  if(fun.typ==7){
    signalfun <- function(x){(2*rbinom(1,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2))}
  }
  if(fun.typ==8){
    signalfun <- function(x){(x > 0.5)}
  }
  signal <- rep(0,n)
  sigmax <- max(signalfun(refx))
  sigmin <- min(signalfun(refx))
  sigmean <- (sigmax+sigmin)/2
  for(i in 1:n){
    if(u[i] > signal.thres){
      signal[i] <- signalfun(x[i]) - sigmean
    } else {
      signal[i] <- 0
    }
  }
  noise <- rep(0,n)
  for(i in 1:n){
    noise[i] <- rnorm(1, sd = noise.level[i] * sd(signal))
  }
  y <- signal + noise
  if(sim.typ==2){
    z <- cbind(z1,z2)
  }
  if(sim.typ==3){
    x <- cbind(x1,x2)
  }
  return(list(x=x,y=y,z=z,signal=signal,noise=noise))
}
