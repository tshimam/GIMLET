gimlet <- function(X, Y, Z, r=1, delta=0.3, nperm=1000, type=c("continuous","binary")){

  type <- match.arg(type)
  if(is.vector(X)) X <- matrix(X,ncol=1)
  if(is.vector(Y)) Y <- matrix(Y,ncol=1)
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)

  if(type=="continuous"){
    result <- gimletCpp(X, Y, Z, r, delta, nperm)
  } else if(type=="binary"){
    result <- gimletCpp2(X, Y, Z, nperm)
  }

  return(result)

}

ldcorAll <- function(X, Y, Z, M, delta=0.3, nperm=1000, type=1){

  if(is.vector(X)) X <- matrix(X,ncol=1)
  if(is.vector(Y)) Y <- matrix(Y,ncol=1)
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)
  if(is.vector(M)) M <- matrix(M,ncol=1)

  result <- ldcorAllCpp(X, Y, Z, M, delta, nperm, type)

  if(nperm >= 1e4){
    each.p.value <- result$each.p.value
    ldcor <- result$ldcor
    null.each.ldcor <- result$null.each.ldcor
    corrected.each.p.value <- each.p.value
    for(i in 1:nrow(ldcor)){
      for(j in 1:ncol(ldcor)){
        if(each.p.value[i,j] >= 1 / (1 + nperm)) next
        tmp.null.each.ldcor <- null.each.ldcor[i,j,]
        thres <- quantile(tmp.null.each.ldcor,0.995)
        upper.null.each.ldcor <- tmp.null.each.ldcor[tmp.null.each.ldcor > thres] - thres
        lambda <- 1 / mean(upper.null.each.ldcor)
        corrected.each.p.value[i,j] <- exp(- lambda * ldcor[i,j]) / exp(- lambda * thres) * 0.005
      }
    }
    result$corrected.each.p.value <- corrected.each.p.value
  }

  return(result)

}
