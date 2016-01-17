#' Genome-wide identification of modulators with local energy test statistics
#'
#' @details This function calculates the local estimator of distance correlation between random vectors X and Y of arbitorary dimension given a specific value of another random vector Z and tests the statistical siginificance of the correlation differences between two different groups of Z.
#' @param X data matrix of X
#' @param Y data matrix of Y
#' @param Z data matrix of Z
#' @param r number of samples in each of two groups with high and low values of Z
#' @param delta tuning parameter controlling the smoothness of the local estimator
#' @param perm number of permutations for calculating p-values
#' @examples 
#' data <- sim(500); result <- gimlet(data$x,data$y,data$z,1,0.3,1000)
#' @return a list consisting of
#' \item{pi}{the approximation of pi}
#' \item{n}{the number of trials}
#' \item{hits}{the number of hits}
#' \item{modulator}{Selected value of Z}
#' \item{ldcor}{Local estimator of distance correlation between X and Y given selected values of Z}
#' \item{global.score}{Test statistic evaluating whether the correlation difference between X and Y depends on Z}
#' \item{each.p.value}{p-value in statistical test under the null hypothesis H_0: the local estimator between X and Y given Z is equal to 0}
#' \item{global.p.value}{p-value in statistical test under the null hypothesis H_0: the local estimator betwee X and Y does not depend on Z}
#' \item{null.global.score}{Test statistics under the null hypothesis}
#' \item{null.each.ldcor}{Local estimators under the null hypothesis}
#' \item{corrected.each.p.value}{Corrected p-value in statistical test under the null hypothesis H_0: the local estimator between X and Y given Z is equal to 0 with the distribution tail approximation}
#' \item{corrected.global.p.value}{Corrected p-value in statistical test under the null hypothesis H_0: the local estimator betwee X and Y does not depend on Z  with the distribution tail approximation}
#' @author Teppei Shimamura & Yusuke Matsui
#' @references Teppei Shimamura, Yusuke Matsui, Taisuke Kajino, Satoshi Ito, Takashi Takahashi and Satoru Miyano. Genome-wide identification of biological modulators with local energy statistics. submitted.
#' @export
gimlet <- function(X, Y, Z, r=1, delta=0.3, perm=1000, type=1){

  if(is.vector(X)) X <- matrix(X,ncol=1)
  if(is.vector(Y)) Y <- matrix(Y,ncol=1)
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)

  result <- gimletCpp(X, Y, Z, r, delta, perm, type)

  if(perm >= 1e4){
    global.score <- result$global.score
    global.p.value <- result$global.p.value
    if(global.p.value <= 1 / (1 + perm)){
      null.global.score <- result$null.global.score
      thres <- quantile(null.global.score,0.995)
      upper.null.global.score <- null.global.score[null.global.score > thres] - thres
      lambda <- 1 / mean(upper.null.global.score)
      corrected.global.p.value <- exp(- lambda * global.score) / exp(- lambda * thres) * 0.005
    } else {
      corrected.global.p.value <- global.p.value
    }

    ldcor <- result$ldcor
    each.p.value <- result$each.p.value
    null.each.ldcor <- result$null.each.ldcor
    corrected.each.p.value <- each.p.value
    for(i in 1:nrow(ldcor)){
      for(j in 1:ncol(ldcor)){
        if(each.p.value[i,j] >= 1 / (1 + perm)) next
        tmp.null.each.ldcor <- null.each.ldcor[i,j,]
        thres <- quantile(tmp.null.each.ldcor,0.995)
        upper.null.each.ldcor <- tmp.null.each.ldcor[tmp.null.each.ldcor > thres] - thres
        lambda <- 1 / mean(upper.null.each.ldcor)
        corrected.each.p.value[i,j] <- exp(- lambda * ldcor[i,j]) / exp(- lambda * thres) * 0.005
      }
    }
    result$corrected.global.p.value <- corrected.global.p.value
    result$corrected.each.p.value <- corrected.each.p.value
  }

  return(result)

}
