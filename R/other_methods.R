Atilde <- function(d) {
  ## d is a distance matrix or distance object
  ## compute u-centered distance matrix
  ## denoted Atilde in Szekely and Rizzo's AS paper (2014)
  d <- as.matrix(d)
  n <- nrow(d)
  if (n != ncol(d)) stop("Argument d should be distance")
  m <- rowSums(d) / (n-2)
  M <- sum(d) / (n-1) / (n-2)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  A <- b + M  #same as plain A
  return(A)
}

pdcor <- function(x, y, z,distance=FALSE) {
  ## compute bias corrected distance correlation
  ## attempt to check if distance flag is valid
  if (distance==FALSE) {
    if (class(x)=="dist" || class(y)=="dist")
      stop("distance==FALSE but argument is a dist object")
    x <- as.matrix(dist(x))
    y <- as.matrix(dist(y))
    z <- as.matrix(dist(z))
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    if (distance == TRUE)
      if (!isSymmetric(x) || !isSymmetric(y))
        stop("distance==TRUE but matrices non-symmetric")
  }
  n <- NROW(x)
  AA <- Atilde(x)
  BB <- Atilde(y)
  CC <- Atilde(z)
  normAA <- sqrt(sum(AA*AA) / n / (n-3))
  normBB <- sqrt(sum(BB*BB) / n / (n-3))
  normCC <- sqrt(sum(CC*CC) / n / (n-3))
  RXY <- sum(AA*BB) / n / (n-3) / normAA / normBB
  RXZ <- sum(AA*CC) / n / (n-3) / normAA / normCC
  RYZ <- sum(BB+CC) / n / (n-3) / normBB / normCC
  r <- (RXY - RXZ*RYZ) / sqrt(1 - RXZ^2) / sqrt(1 - RYZ^2)
  return(r)
}

pdcov.test <- function(x, y, z,distance=FALSE,B=999) {
  ## compute bias corrected distance correlation
  ## attempt to check if distance flag is valid
  if (distance==FALSE) {
    if (class(x)=="dist" || class(y)=="dist")
      stop("distance==FALSE but argument is a dist object")
    x <- as.matrix(dist(x))
    y <- as.matrix(dist(y))
    z <- as.matrix(dist(z))
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    z <- as.matrix(z)
    if (distance == TRUE)
      if (!isSymmetric(x) || !isSymmetric(y))
        stop("distance==TRUE but matrices non-symmetric")
  }
  n <- NROW(x)
  AA <- Atilde(x)
  BB <- Atilde(y)
  CC <- Atilde(z)
  U = Atilde(AA - sum(AA*CC)/sum(CC*CC) * CC)
  V = Atilde(BB - sum(BB*CC)/sum(CC*CC) * CC)
  obj <- dcov.test(U,V,R=B)
  return(obj)
}

require(infotheo)

mindy <- function(X,Y,Z,B=1000,method="emp",percent=0.35){
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)
  X <- discretize(X)[,1]
  Y <- discretize(Y)[,1]
  n <- length(X)
  m1 <- round(n*percent,0)
  m2 <- round(n*(1-percent),0)
  id1 <- sort.list(Z)[1:m1]
  id2 <- sort.list(Z)[m2:n]
  score <- abs(mutinformation(X[id2],Y[id2],method)-mutinformation(X[id1],Y[id1],method))
  null.score <- rep(0,B)
  for(i in 1:B){
    id1 <- sample(1:n,m1,replace=FALSE)
    id2 <- sample(1:n,m2,replace=FALSE)
    score0 <- abs(mutinformation(X[id2],Y[id2],method)-mutinformation(X[id1],Y[id1],method))
    null.score[i] <- score0
  }
  p.value <-  (sum(score<null.score) + 1) / (B + 1)
  return(list(score=score,perm.score=null.score,p.value=p.value))
}

gem <- function(X,Y,Z,percent=0.3){
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)
  m <- length(X)
  if(ncol(Z)>=2) Z <- sqrt(apply(Z^2,1,sum))
  X2 <- Y2 <- Z2 <- rep(NA,m)
  X2[X>=quantile(X,probs=(1-percent))] <- 1
  X2[X<=quantile(X,probs=percent)] <- 0
  Y2[Y>=quantile(Y,probs=(1-percent))] <- 1
  Y2[Y<=quantile(Y,probs=percent)] <- 0
  Z2[Z>=quantile(Z,probs=(1-percent))] <- 1
  Z2[Z<=quantile(Z,probs=percent)] <- 0
  f <- array(0,c(2,2,2))
  for(i in 1:2) for(j in 1:2) for(k in 1:2) f[i,j,k] <- sum(Z2==(i-1)&X2==(j-1)&Y2==(k-1),na.rm=TRUE)
  n <- apply(f,c(1,2),sum)
  p <- f[,,2] / n
  if(any(n==0)) p[n==0] <- 0
  alpha.c <- p[1,1]
  alpha.f <- p[1,2] - p[1,1]
  alpha.m <- p[2,1] - p[1,1]
  beta.f <- p[2,2] - p[2,1]
  beta.m <- p[2,2] - p[1,2]
  gamma <- p[2,2] - p[1,2] - p[2,1] + p[1,1]
  p.alpha.f <- (f[1,2,2] + f[1,1,2]) / (n[1,2] + n[1,1])
  p.alpha.m <- (f[2,1,2] + f[1,1,2]) / (n[2,1] + n[1,1])
  p.beta.f <- (f[2,2,2] + f[2,1,2]) / (n[2,2] + n[2,1])
  p.beta.m <- (f[2,2,2] + f[1,2,2]) / (n[2,2] + n[1,2])
  p.gamma <- (f[1,1,2] + f[1,2,2] + f[2,1,2] + f[2,2,2]) / (n[1,1] + n[1,2] + n[2,1] + n[2,2])
  var.alpha.f <- p.alpha.f * (1 - p.alpha.f) * (1 / n[1,2] + 1 / n[1,1])
  var.alpha.m <- p.alpha.m * (1 - p.alpha.m) * (1 / n[2,1] + 1 / n[1,1])
  var.beta.f <- p.beta.f * (1 - p.beta.f) * (1 / n[2,2] + 1 / n[2,1])
  var.beta.m <- p.beta.m * (1 - p.beta.m) * (1 / n[2,2] + 1 / n[1,2])
  var.gamma <- p.gamma * (1 - p.gamma) * (1 / n[1,1] + 1 / n[1,2] + 1 / n[2,1] + 1 / n[2,2])
  # z <- function(x){exp(-x^2)}
  # 1 - 2 / sqrt(pi) * integrate(z, lower=0,upper= abs(beta.m / sqrt(2*var.beta.m)))$value
  # is equivalent to
  p.value.alpha.f <- pnorm(abs(alpha.f),0,sqrt(var.alpha.f),lower.tail=FALSE)*2
  p.value.alpha.m <- pnorm(abs(alpha.m),0,sqrt(var.alpha.m),lower.tail=FALSE)*2
  p.value.beta.f <- pnorm(abs(beta.f),0,sqrt(var.beta.f),lower.tail=FALSE)*2
  p.value.beta.m <- pnorm(abs(beta.m),0,sqrt(var.beta.m),lower.tail=FALSE)*2
  p.value.gamma <- pnorm(abs(gamma),0,sqrt(var.gamma),lower.tail=FALSE)*2
  flag <- FALSE
  if(abs(beta.m)>abs(alpha.m)|sign(beta.m)== - sign(alpha.m)) flag <- TRUE
  result <- list(alpha.c=alpha.c,alpha.f=alpha.f,alpha.m=alpha.m,
                       beta.f=beta.f,beta.m=beta.m,gamma=gamma,
                       p.value1=p.value.beta.m,p.value2=p.value.gamma,flag=flag)
  return(result)
}
